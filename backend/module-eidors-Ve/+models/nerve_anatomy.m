
function nerve_anatomy (varargin)
% nerve_anatomy(geometry, ...) runs the electroanatomical model implemented
% using EIDORS for electrical stimulation & recording from a peripheral 
% nerve, including meshing the model in gmesh and gnerating any neccessary
% thin-layers
% 
% This is iteration #3 of this code, which has fixed issues with lost
% mesh points and uses ./+mesh/pelvic_nerve.geo.template
% (future versions will have more ways of customising the models). 
% 
% "geometry" is a struct which defines the stimulating / recording
% electrode geometry, some meshing parameters, and which fascicle set to
% use (from the default splines.dat, see read_dat_file). If "geometry" is
% not supplied, it's read from tools.cache('PATH','elec-geom.mat')
% 
% Some additional -options can be specified: 
% -regenerate: if set, the mesh is regenerated even if
%              tools.cache('path')/pelvic_nerve-thin.msh.mat exists. 
% -mesh:     generate mesh only and exit without running EIDORS. 
% -fascicle: which fascicles from file get simulated (default: 1-5) 
% -stimulus: compute the electrode-tissue relationship for electrical 
%            stimuli (if not set, this code computes recording sensitivity)
%     -mono: if set, stimuli are monopolar (default: sequential bipoles)
% -reference: default [ COUNTER = 1mm, -1, d=2   GROUND = 1mm, -1, d=3 ] 
% -common: use the COUNTER for both COUNTER and GROUND
% -output: sets output file, default '~/eidors/sensitivity (%d).mat'
% 
% Calvin Eiber 25-Apr-2020

%% Parse model inputs

tools.setupEIDORS;

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

% save file to be read by insert_gmsh_electrodes during compilation:
% @[mesh.insert_gmsh_electrodes('-import','elec_geom.mat')] 
% if exist('e','var'), save elec-geom.mat -struct e , end %#ok<USENS>

e = get_input_structure(varargin{:}); 

if any(named('-out')), output_name = get_('-out');
elseif isfield(e,'output'), output_name = e.output;
else output_name = ''; 
end

% if output_name exists but does not end in .mat: 
if ~isempty(output_name) && isempty(regexp(output_name,'\.mat$','once'))
    output_name = sprintf('%s (%s).mat', tools.file('eidors~/sensitivity'), ...
                                         output_name);
elseif isempty(output_name) 
    output_name = tools.file('get','eidors~/sensitivity (%d)', 'next');
end

if ~exist(fileparts(output_name),'dir')
    warning('ViNERS:nerveAnat:missingOutputDir', ...
            'requested output dir %s did not exist. This often %s (%s).', ...
            fileparts(output_name), ...
            'indicates a problem with the current subject in tools.file', ...
            tools.file('T',tools.file('sub~')))
    mkdir(fileparts(output_name));
end

%% Generate requested mesh
m = generate_nerve_mesh(get_,named, e); 
m = configure_array_electrodes(m); % Configure electrodes. 

%% Show the volume and electrode pattern
preview_eidors_mesh(m)

%% Preview check! [fast-ish]
if any(named('-mesh')), return, end

%% Construct virtual point electrodes and stimulation pattern for computation of sensitivity

disp('Constructing current pattern ... ')
% note: this code depends on exact m.object_names 

m.electrode(cellfun(@isempty,{m.electrode.name})) = [];
nE = numel(m.electrode);

clf, show_fem(m), hold on
h = get(gca,'Children'); 
set(h,'EdgeAlpha',0.3)
set(h,'FaceAlpha',0.1)

obj_idx_ = @(n) unique( m.elems(cat(1,m.object_id{contains(m.object_name,n)}),:));

% e_xyz = mean(m.nodes(m.electrode(floor(end/2)).nodes,:)); 
% e_ref = m.electrode(floor(end/2)).name;
% [~,idx] = sort( sqrt(sum((m.nodes-e_xyz).^2,2)) ); 

if any(named('-fasc'))
     fascicle_list = get_('-fasc');
else fascicle_list = 1:sum(strncmp(m.object_name,'Fascicle',8)); 
end

for f_id = fascicle_list % Each fascicle, broken into chunks. 
  
    chunk_size = 1000; 
    fascicle_id = sprintf('Fascicle%d', f_id);

    [~,idx] = sort( m.nodes(:,1) ); % all points, left-to-right
    i_sim = obj_idx_(fascicle_id);  % FascicleN and P_FascicleN    
    
    idx(~ismember(idx,i_sim)) = []; % idx becomes i_sim but sorted
    
    % idx(m.nodes(idx,1) < 0) = []; Array not garunteed to be symmetrical
    idx(m.nodes(idx,2) < 0) = [];   % remove negative Y
    idx(ismember(idx,m.boundary(:))) = []; % remove boundary
    
    idx_full = idx; 
    
    if numel(idx) < 1800, chunk_size = 1800; 
    else n_chunk = round(numel(idx)/chunk_size);
         chunk_size = ceil(numel(idx)/n_chunk);
    end
    
    %% CORE loop (only relevent if -direct set)
    for chunk = 1:max(1,ceil(numel(idx_full) / chunk_size))
        %% Do in lots of 1000 points, keeps the slowness managable.
        m.electrode(cellfun(@isempty,{m.electrode.name})) = [];

        idx = (1:chunk_size) + chunk_size * (chunk-1);
        idx(idx > numel(idx_full)) = []; 
        idx = idx_full(idx); 

        if isempty(idx), continue, end
        
        if any(named('-set-zc'))
          [m.electrode.z_contact] = deal(get_('-set-zc'));
        end
        
        nP = numel(idx);
        ref = nE + nP + 1; 

        v_ = @(x) reshape(x,1,[]);
        meas_pattern = sparse(v_((1:nE)' * [1 1]), [1:nE ref*ones(1,nE)], ...
                                                [ones(1,nE) -ones(1,nE)]);
        m.stimulation = struct([]);
        
        if any(named('-stim')) || ~any(named('-direct'))
          ref = nE + 1;
          meas_pattern = sparse(v_((1:nE)' * [1 1]), [1:nE ref*ones(1,nE)], ...
                                                  [ones(1,nE) -ones(1,nE)]); % Updated reference
          for ii = 1:nE 
              m.stimulation(ii).stimulation = 'Amp';
              if any(named('-mono')) || nE == 1 || ~any(named('-stim'))
                   stim_pattern = sparse([ii ref],[1 1],[1 -1],nE+1,1);
              else stim_pattern = sparse([ii mod(ii,nE)+1],[1 1],[1 -1],nE+1,1);
              end
              m.stimulation(ii).stim_pattern = stim_pattern;
              m.stimulation(ii).meas_pattern = meas_pattern(stim_pattern(1:end-1) == 0,:);
          end        
        elseif any(named('-veri')), m = mk_verification_pattern(m);
        else 
          m.stimulation(nP).stimulation = 'Amp';
          m.stimulation(nP).stim_pattern = sparse([nP+nE ref],[1 1],[1 -1]);
          m.stimulation(nP).meas_pattern = meas_pattern; 
          for ii = 1:nP % Insert point current sources 
            m.electrode(ii+nE).nodes = idx(ii);
            m.electrode(ii+nE).z_contact(ii+1) = 1;
            m.stimulation(ii).stimulation = 'Amp';
            m.stimulation(ii).stim_pattern = sparse([ii+nE ref],[1 1],[1 -1]);
            m.stimulation(ii).meas_pattern = meas_pattern;
          end
          plot3(m.nodes(idx(1:nP),1), m.nodes(idx(1:nP),2), m.nodes(idx(1:nP),3), '.')
          pause(0.1)
        end
        
        %% Set up reference electrode
        m = setup_ref_electrodes(m,ref,chunk);

        %% Set conductivities and run model
        
        fprintf('Running field simulation [%s chunk %d] ... \n', fascicle_id, chunk)
        if any(named('-sigma')) % Import volume conductivities from fitted transimpedance data
             em = build_forward_model(m,get_('-sigma')); % use custom anisotropic conductivities
        else em = build_forward_model(m); % set up non-homogenous anisotropic conductivity
        end 
        
        if any(named('-eit-model')) % Import volume conductivities from fitted transimpedance data
            em = import_fitted_conductivities(em, get_('-eit-model')); % set up non-homogenous anisotropic conductivity
        end
        
        if exist('elem_data_isotropic','var')
            em.elem_data = elem_data_isotropic;
        end
        
        if any(named('-stimulus')) || ~any(named('-direct'))

           em.fwd_solve.get_all_meas = 1; % Get ALL points
           v_meas = fwd_solve(em); % solve EIDORS electric field simulation 
           v_meas = v_meas.volt;             

           ok = cat(1,m.object_id{strncmp(m.object_name,'P_Fa',3)});
           if ~any(ok)
             ok = cat(1,m.object_id{strncmp(m.object_name,'Fascicle',3)});
           end

           if any(ok)
             ok = unique(m.elems(ok,:));
             scatter3(m.nodes(ok,1), m.nodes(ok,2), m.nodes(ok,3), 10, ...
                              v_meas(ok,1),'o'), colormap(tools.magma)
             caxis(quantile(v_meas(ok,1),[.01 .99]));
           end
           
           if isempty(output_name), output_name = tools.file('get', ... 
                                'sub~/eidors/stimulus (%d).mat','next'); 
           elseif contains(output_name,'sensitivity')
             output_name = strrep(output_name,'sensitivity','stimulus');
           elseif exist(output_name,'file')
             output_name = tools.file('get',output_name,'next'); 
           end
      
           m.name = 'Extracellular potential simulation';
           info = e;
           model = m;
           v_extracellular = v_meas;
           save(output_name,'info','model','v_extracellular')
           
           if ~any(named('-direct')) && ~any(named('-stimulus'))           
             % see https://en.wikipedia.org/wiki/Reciprocity_(electromagnetism)#Reciprocity_for_electrical_networks
               result.info = e;
               result.model = m;
               result.v_extracellular = v_meas;
               result = convert_stim2fascicles(result); 
               result = rmfield(result,'v_extracellular');
                 
               output_name = strrep(output_name,'stimulus','sensitivity');
               save(output_name,'-struct','result');     
           end
           
           return
        else
            
          v_meas = fwd_solve(em); % solve EIDORS electric field simulation
          v_meas = reshape(v_meas.meas,nE,[]);        
          disp('Saving to "eidors-im_pot.mat" ... ')
          save(tools.file('get',tools.cache('PATH', ... 
                          'eidors-im_pot_%s_chunk%02d (%%d).mat'), ... 
                          'next',fascicle_id,chunk),'m','v_meas')
        end
        
        if any(isnan(v_meas(:)))
          eidors_cache clear
          error('Something went awry, %0.1f%% of v_meas came out NaN', ...
                                       100*mean(isnan(v_meas(:))))
        end
        if any(named('-veri')), break, end % for each fascicle
    end % chunk loop

    % % Debug view using EIDORS plot function: 
    % clf, show_slices(img_bkgnd,[0 inf inf; inf 0 inf; inf inf 0])
end % fascicle_ID loop

%%
fprintf('Done!\n')

if any(named('-fix-m')), coalesce_EIDORS(old_eidors_file,output_name);  
else coalesce_EIDORS(output_name,'new');
end

return % everything after this is visualisation

%%
figure(1) %#ok<UNRCH>
clf, show_fem(m), hold on

nE = find(cellfun(@isempty,{m.electrode.name}),1)-1;
nP = size(v_meas,2);
ee = 2; 

h = get(gca,'Children'); 
hp = arrayfun(@(o) isa(o,'matlab.graphics.primitive.Patch'),h); 
delete(h(~hp))
h = h(hp);

set(h,'EdgeAlpha',0.3)
set(h,'FaceAlpha',0.5)
set(h(1),'FaceAlpha',0)

meas_xyz = m.nodes([m.electrode(nE+(1:nP)).nodes],:);

% patch('Faces',model.elems(carrier,:),'Vertices',model.nodes,'EdgeColor','r','FaceAlpha',0)

scatter3(meas_xyz(:,1),meas_xyz(:,2),meas_xyz(:,3), [], v_meas(ee,:),'s','filled')
title(m.electrode(ee).name)

caxis([min(v_meas(ee,:)) max(v_meas(ee,:))])
colormap(tools.magma), colorbar('SouthOutside')

zoom = 0.02 * [1 -1; -1 1] + [1 0; 0 1];
axis(axis * blkdiag(zoom,zoom,zoom)), grid on
set(gca,'Color','none')

clear h hp zoom ee

% axis([0 2.2 -1 1 -2 2])

%% Distance - sensitivity measure 

figure(2), clf, hold on

for ee = 1:nE 
    e_xyz = mean(m.nodes(m.electrode(ee).nodes,:),1);
    e_dist = sqrt(sum((e_xyz - meas_xyz).^2,2));
    plot(e_dist,v_meas(ee,:),'o','MarkerSize',4)
end

tools.tidyPlot
xlabel('Distance from electrode, mm')
ylabel('Sensitivy, uV / uA')
clear ee d0 

% set(gca,'YScale','log')
% The response (uV/uA) seems to be linear w.r.t. conductivity as expected
return

%% Load the array configuration from .mat or .json file; check defaults 
function e = get_input_structure(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

is_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if nargin > 0 && isstruct(varargin{1}), file = varargin{1}; 
elseif nargin > 0 && exist(varargin{1},'file'), file = varargin{1}; 
elseif any(named('-geom')), file = get_('-geom'); 
elseif exist(tools.cache('path','elec-geom.mat'),'file')
        file = tools.cache('path','elec-geom.mat'); 
elseif exist(tools.cache('path','elec-geom.json'),'file')
        file = tools.cache('path','elec-geom.json');
else    file = models.electrode_array(varargin{:});
end

if isstruct(file), e = file; 
else   
  if is_ext_(file,'.mat'),      e = load(file); 
  elseif is_ext_(file,'.json'), e = tools.parse_json(file);
  elseif is_ext_(file,'.xml'),  e = tools.parse_xml(file); 
  else error('unknown filetype on file "%s", expected {.mat, .json, .xml}', file)
  end

  if ~any(named('-q')), fprintf('Loading %s\n', file); end
  % CE 9-4-21: I haven't tested the .json or .xml here yet.   
end

if isfield(e,'mesh') % move these fileds to 'top-level'
  for fi = fieldnames(e.mesh)'
    e.(fi{1}) = e.mesh.(fi{1}); 
  end
end

%% Add some default fields 

if ~isfield(e,'MeshLengthMax')
    e.MeshLengthMax = 0.1;   
    warning('ViNERS:nerveAnat:default', ...
            'using default MeshLengthMax, %g mm', e.MeshLengthMax)
end
% if ~isfield(e,'SplineIndex'),    e.SplineIndex = 100;     end
if ~isfield(e,'Perineurium_mm') 
  if isfield(e,'nerve') && isfield(e.nerve,'Perineurium_mm')
       e.Perineurium_mm = e.nerve.Perineurium_mm;
  elseif isfield(e,'nerve') && isfield(e.nerve,'Perineurium_um')
       e.Perineurium_mm = e.nerve.Perineurium_um / 1e3;
  else e.Perineurium_mm = 1e-3; 
       warning('ViNERS:nerveAnat:default', ...
               'using default perineurium, %g um', e.Perineurium_mm*1e3)
  end
end

if ~any(named('-no-cache'))
    save(tools.cache('path','elec-geom.mat'),'-struct','e')
end


return

%% Collect intermediate results into final output structure
function out = coalesce_EIDORS(outname,inputs)

default_in = '*(1).mat'; % <<< Which EIDORS sim do we want? 
cache_path = tools.cache('path'); % <<< Where do we find it?

if nargin < 2, inputs = default_in; end
if nargin < 1 || isempty(outname)
  outname = tools.file('get','out~/sensitivity (%d).mat','next');
end

if isstruct(outname) % append mode   
  out = outname;
  outname = inputs;
  inputs = 'new';
else out = []; 
end

if any(inputs == '*'), list = dir([cache_path inputs]); 
elseif strncmpi(inputs,'new',3)    
    list = dir([cache_path '*.mat']); % everything
    [~,ff] = max([list.datenum]); % newest
    
    inputs = regexprep(list(ff).name,'[^\-\(]+\d+','*'); % make filename less unique
    list = dir([cache_path inputs]); % get matching simpler filename
end

if strcmp(inputs,'new;GETNAME'), out = outname; return, end

%%
save_stimpattern = 0; 

for ff = 1:length(list)
    
    f_id = regexp(list(ff).name,'Fascicle\d+','match','once');
    
    if isempty(f_id), 
        warning('ViNERS:EidorsParse', ...
                'Unclear how to parse "%s"',list(ff).name)
    end
    
    % Load chunk
    in = load([cache_path list(ff).name]);    
    idx = cellfun(@length, {in.m.electrode.nodes}) == 1; %% Point sources
    xyz = [in.m.electrode(idx).nodes];
   
    if isempty(out), % Initialise
        out = struct;
        out.model = in.m; 
        out.model.electrode(idx) = []; 
        out.model.stimulation = []; 
        
        if exist('elec-geom.mat','file')
          out.info = load('elec-geom.mat'); 
        elseif exist(tools.cache('path','elec-geom.mat'),'file')
          out.info = load(tools.cache('path','elec-geom.mat'));
        else
          out.info.SplineIndex = 100; % default fascicle
        end
        
        if ~isempty(which('mk_inhomogenous')) % on path?
            out.info.conductivity = build_forward_model;
        end
        
        ret = cellfun(@(x) find(x==-1), {in.m.stimulation.stim_pattern});
        save_stimpattern = save_stimpattern | numel(unique(ret)) > 2; 
    end
    if ~isfield(out,f_id)
       out.(f_id).pot = [];
       out.(f_id).idx = []; 
       if save_stimpattern, out.(f_id).stim = {}; end
    end        
    out.(f_id).pot = [out.(f_id).pot in.v_meas]; 
    out.(f_id).idx = [out.(f_id).idx xyz]; 
    
    if save_stimpattern
       out.(f_id).stim = [out.(f_id).stim {in.m.stimulation}]; 
    end
end

out.utils.x_ = '@(i) out.model.nodes(i,1)';
out.utils.y_ = '@(i) out.model.nodes(i,2)';
out.utils.z_ = '@(i) out.model.nodes(i,3)';
out.utils.fac_ = '@(i) out.(sprintf(''Fascicle%d'',i)).idx'; 

if ~exist(fileparts(outname),'dir'), mkdir(fileparts(outname)), end
disp(['Saving ' outname])
save(outname,'-struct','out')

return
%% Quick test

disp('Running quick debug plot') %#ok<UNRCH>

for f = fieldnames(out.utils)' % Create utility local functions
    eval(sprintf('%s = %s;',f{1},out.utils.(f{1})));
end

nF = nanmax(str2double(regexp(fieldnames(out),'\d+','match','once')));

for ff = 1:nF
    if ~isfield(out,sprintf('Fascicle%d',ff)), continue, end    
    plot(z_(fac_(ff)),y_(fac_(ff)),'.')
    hold on, axis equal
end

tools.tidyPlot

%% Generate mesh from specified inputs 
function m = generate_nerve_mesh(get_,named, e)

PN_mesh = tools.cache('PATH','nerve+array'); 
PN_mesh = @(ext)[PN_mesh ext]; % output file

if any(named('-regenerate')) || ~exist(PN_mesh('-thin.msh.mat'),'file')

  if any(named('-template')), geo_template = get_('-template');
  elseif isfield(e,'array') && isfield(e.array,'Template')
       geo_template = e.array.Template;
  elseif isfield(e,'array') && isfield(e.array,'template')
       geo_template = e.array.template;
  else % Not supplied as input argument, go find it      
    geo_template = tools.file('get','mesh~\*.geo.template','-q');
    if isempty(geo_template)
        geo_template = 'source~\array\planar-array.geo.template'; 
    end
  end
  
  if any(geo_template == '~'), geo_template = tools.file(geo_template); end
  if ~exist(geo_template,'file')
      geo_template = [geo_template '.geo.template'];
  end
  t = tic;
  
  make_opts = {'-usev','-output',PN_mesh('.geo'), ...
                  'virtual_thinlayer','exterior_len',e.MeshLengthMax};
  if isfield(e,'array') && isfield(e.array,'carrier') % case-sensitive
      more_opts = tools.opts_to_args(e.array,'carrier','--s2a-keep');
      make_opts = [make_opts more_opts];
  end
  
  if any(named('-make-o')), more_opts = get_('-make-o');
    if ~iscell(more_opts), more_opts = {more_opts}; end
    make_opts = [make_opts more_opts];
  end
  
  if ~exist(PN_mesh('.geo'),'file')
    tools.makeFromTemplate(geo_template, make_opts{:})
  else   warning('ViNERS:existing_geo_file', ...
                 'Using existing .geo file (%s). %s', PN_mesh('.geo'), ...
                 'If this is not intended, call tools.cache(''reset'')')
  end
  
  % use gmsh.exe to convert the .geo file to a .msh file:
  if ~exist(PN_mesh('.msh'),'file')  
    path_to_gmsh = tools.configuration('gmsh'); % 'C:\Program Files\gmesh\gmsh.exe';
    system(sprintf('"%s" "%s" -format msh41 -3',path_to_gmsh,PN_mesh('.geo'))); 
  end

  if ~exist(PN_mesh('.msh'),'file')
      
    if isunix, system(['gmsh ' PN_mesh('.geo')])
    else               winopen(PN_mesh('.geo'))
    end
    error('Pipeline crashed during mesh generation')
  end
  
  % Add thinlayer and convert to EIDORS forward model (takes a minute or two)
  
  m = mesh.make_gmsh_thinLayer(PN_mesh('.geo'), e.Perineurium_mm );
  fprintf('Elapsed time is %.0f:%02.3f\n', floor(toc(t)/60), mod(toc(t),60))  
else
  warning('ViNERS:existing_mat_file', ...
                 'Using existing mat file (%s). %s', PN_mesh('-thin.msh.mat'), ...
                 'If this is not intended, call tools.cache(''reset'')')
  
  m = load(PN_mesh('-thin.msh.mat'));
end
  
m.boundary_numbers = ones(size(m.boundary(:,1)));

%% Add internal electrodes to model 
function m = configure_array_electrodes(m)
% configure_electrodes introduces the electrodes to the model as empty
% spaces (voids) - needed for internal electrodes with geometry in EIDORS

if ~isfield(m,'electrode') % Find electrode in model and 
                           % eliminate corresponding elements 
    % Electrodes already exist as voids, dig them out here. 
    error('Locate electrode voids code removed from an older version')
end

for ee = 1:numel(m.electrode)
    sel = all(ismember(m.elems',m.electrode(ee).nodes))';    
    m.elems(sel,:) = 0;
end

rzc = cumsum(m.elems(:,1) == 0); % running zero count

for ii = 1:numel(m.object_id) % Re-index object map
    oid = m.object_id{ii}; 
    sel = m.elems(oid,1) ~= 0;
    m.object_id{ii} = oid(sel) - rzc(oid(sel));
end

m.elems(m.elems(:,1) == 0,:) = [];     
m.boundary = find_boundary(m);
m.boundary_numbers = ones(size(m.boundary(:,1)));

for ee = 1:numel(m.electrode) % Get boundary numbers

    sel = all(ismember(m.boundary',m.electrode(ee).nodes));
    m.boundary_numbers(sel) = ee + 1; 
end
clear ee ii sel oid rcz

%% Clean up artifacts from electrode introduction process

nE = numel(m.electrode); 
for ee = 1:nE
    
    ok = ismember(m.electrode(ee).nodes,m.boundary);
    if all(ok), continue, end
    fprintf('[%c%d nodes removed from %s]%c\n',8,sum(~ok),m.electrode(ee).name,8)
    bad = m.electrode(ee).nodes(~ok);
    m.electrode(ee).nodes(~ok) = []; 
    
    % re-index everything again?    
    e_bad = any(ismember(m.elems,bad),2); 
    if ~any(e_bad)
      reindex = arrayfun(@(x) x-sum(bad<=x), 1:max(m.elems(:)));  
      reindex(bad) = -1; 
      m.nodes(bad,:) = [];
      m.elems = reindex(m.elems);
      m.boundary = reindex(m.boundary);
     
      for ii = 1:numel(m.electrode)
        m.electrode(ii).nodes = reindex(m.electrode(ii).nodes);
      end
%     for ii = 1:numel(m.object_id) % index over elems doesn't change
%       m.object_id{ii} = reindex(m.object_id{ii});
%     end      
    end    
end

%% Need to check for orphans, otherwise sys_mat is singular
indeg = 0*m.nodes(:,1); 
for ii = 1:size(m.elems,1), % Much faster than equivalent arrayfun
  indeg(m.elems(ii,:)) = indeg(m.elems(ii,:)) + 1;
end

if any(indeg == 0)
  
  disp('Cleaning up system matrix ... ')
  
  bad = find(indeg == 0);
  reindex = arrayfun(@(x) x-sum(bad<=x), 1:max(m.elems(:)));  
  m.nodes(bad,:) = [];
  m.elems = reindex(m.elems);
  m.boundary = reindex(m.boundary);
      
  for ii = 1:numel(m.electrode)
    m.electrode(ii).nodes = reindex(m.electrode(ii).nodes);
  end
end

%% Verification pattern
function m = mk_verification_pattern(m) % IN-CONTEXT function

nE = evalin('caller','nE');
fascicle = evalin('caller','idx_full');

%% Simulate semi-random pairs of points in addition to monopolar points. 
%    FOR EACH electrode, 50% are distributred across nodes near other
%    electrodes and 50% are distributed 

nV = min(3*nE,20);
e_nIndex = cell(nE,1); 
v_ = @(x) reshape(x,1,[]);

for ee = 1:nE % Get "close to electrode" set 

  e_xyz = mean(m.nodes(m.electrode(ee).nodes,:),1);
  [~,idx] = sort(sum((m.nodes(fascicle,:) - e_xyz).^2,2));
  e_nIndex{ee} = fascicle(idx(1:min(2*nV,end)));
end

o_nIndex = setdiff(fascicle,cat(1,e_nIndex{:}));
[~,seq] = sort(sum(m.nodes(o_nIndex,:).^2,2));
seq = seq(unique(round(linspace(1,numel(seq),nV*nE))));
o_nIndex = o_nIndex(seq);

E_MAP = unique([cat(1,e_nIndex{:}); o_nIndex]);
for ii = 1:numel(E_MAP) % Insert point current sources 
  m.electrode(ii+nE).nodes = E_MAP(ii);
  m.electrode(ii+nE).z_contact = 1;
end
E_DONE = false(size(E_MAP));

nP = numel(E_MAP);
ref = nE + nP + 1; % Revised reference index
meas_pattern = sparse(v_((1:nE)' * [1 1]), [1:nE ref*ones(1,nE)], ...
                                        [ones(1,nE) -ones(1,nE)]);

assignin('caller','ref',ref);
assignin('caller','nP',nP);
assignin('caller','meas_pattern',meas_pattern);

m.name = 'Verification Test Pattern';
m.stimulation = struct;
m.stimulation(nE).stimulation = 'Amp';
m.stimulation(nE).stim_pattern = sparse([nP+nE ref],[1 1],[1 -1]);
m.stimulation(nE).meas_pattern = meas_pattern; 

o_snk = 1; 
e_snk = 0; 
s = 0; 

while ~all(E_DONE)
  
  e_snk = mod(e_snk,nE)+1; % desynch to make e1-e2, e1-e3, etc. 
  
  any_new = false;
  
  for e_src = 1:nE    
    %% Add elec-elec point pair
    p1 = e_nIndex{e_src}; 
    p2 = e_nIndex{e_snk};
    
    p_done = arrayfun(@(p) E_DONE(E_MAP == p), p1);    
    p1 = p1(find(~p_done,1));
    if isempty(p1), continue, end
    
    p2(p2 == p1) = [];     
    p2 = p2(ceil(rand * numel(p2)));        
    
    s = s+1;     
    s_idx = [find(E_MAP == p1) find(E_MAP == p2) ref];
    m.stimulation(s).stimulation = 'Amp';
    m.stimulation(s).stim_pattern = sparse(s_idx,[1 1 1],[1 -1 0]);
    m.stimulation(s).meas_pattern = meas_pattern;
    
    e_snk = mod(e_snk,nE)+1;
    E_DONE(p1 == E_MAP) = 1;
    any_new = true;
    
    %% Add elec-other point pair
    p1 = e_nIndex{e_src}; 
    p2 = o_nIndex(o_snk);
    
    p_done = arrayfun(@(p) E_DONE(E_MAP == p), p1);    
    p1 = p1(find(~p_done,1));
    if isempty(p1), continue, end
    
    s = s+1;     
    s_idx = [find(E_MAP == p1) find(E_MAP == p2) ref];
    m.stimulation(s).stimulation = 'Amp';
    m.stimulation(s).stim_pattern = sparse(s_idx,[1 1 1],[1 -1 0]);
    m.stimulation(s).meas_pattern = meas_pattern;
    
    o_snk = mod(o_snk,numel(o_nIndex))+1;
    E_DONE(p1 == E_MAP) = 1; 
    E_DONE(p2 == E_MAP) = 1;
    any_new = true;
    
  end
  
  if ~any_new, break, end
end

%%
for ii = 1:nP % Insert monopolar-source simulations 
   
  s = s+1;     
  m.stimulation(s).stimulation = 'Amp';
  m.stimulation(s).stim_pattern = sparse([ii+nE ref],[1 1],[1 -1]);
  m.stimulation(s).meas_pattern = meas_pattern;
end
plot3(m.nodes(E_MAP(1:nP),1), m.nodes(E_MAP(1:nP),2), m.nodes(E_MAP(1:nP),3), '.')
pause(0.1)

%% Add one or two reference electrodes to outer surfaces 
function m = setup_ref_electrodes(m,ref,chunk) % IN-CONTEXT function

named = evalin('caller','named');
get_ = evalin('caller','get_');
nE = evalin('caller','nE');
nP = evalin('caller','nP');

if any(named('-reference')),  ref_electrode = get_('-reference');
elseif any(named('-common')), ref_electrode = [2 -1 2]; 
else ref_electrode = [2 -1 2; 1 -1 3]; 
end

if any(named('-stimulus')), ref_electrode = ref_electrode(end,:); end

d = ref_electrode(1,3); % which side to put on? [xyz] = 1-3
idx = unique(m.boundary(m.boundary_numbers == 1,:));
if ref_electrode(1,2) < 0
     idx(m.nodes(idx,d) > min(m.nodes(idx,d))+10*eps) = []; 
else idx(m.nodes(idx,d) < max(m.nodes(idx,d))-10*eps) = []; 
end

d = ~ismember(1:3,d); 
e_dist = sqrt(sum(m.nodes(idx,d).^2,2)); 
if ~any(e_dist > ref_electrode(1,1))
  warning('ViNERS:pointRefElec','COUNTER (ref 1) is a point electrode')
  [~,idx] = min(e_dist);
else idx(e_dist > ref_electrode(1,1)) = []; 
end

m.electrode(ref).nodes = idx;
m.electrode(ref).z_contact = 1; % changing to 0 results in NaN output
m.gnd_node = m.electrode(ref).nodes(1); 

if chunk == 1 && any(named('-reference'))
  plot3(m.nodes(idx,1),m.nodes(idx,2),m.nodes(idx,3),'k.','MarkerSize',4)
  text(mean(m.nodes(idx,1)),mean(m.nodes(idx,2)),mean(m.nodes(idx,3)),'COUNTER')
end

if any(named('-bipolar'))
  error TODO_construct_bipolar_meas_pattern
elseif size(ref_electrode,1) > 1 
  %% Use a second virtual distant electrode for recording reference
  % The mesh-free simplified model suggusted that the "pedestal"
  % which sensitivity peak was sitting on was be artifactual. I
  % confirmed that using a dual-referece model which I am now
  % making the default. 

  d = ref_electrode(2,3); % which side to put on? [xyz] = 1-3
  idx = unique(m.boundary(m.boundary_numbers == 1,:));
  if ref_electrode(2,2) < 0
       idx(m.nodes(idx,d) > min(m.nodes(idx,d))+10*eps) = []; 
  else idx(m.nodes(idx,d) < max(m.nodes(idx,d))-10*eps) = []; 
  end

  d = ~ismember(1:3,d);
  e_dist = sqrt(sum(m.nodes(idx,d).^2,2)); 
  if ~any(e_dist > ref_electrode(2,1))
    warning('ViNERS:pointRefElec','GROUND (ref 2) is a point electrode')
    [~,idx] = min(e_dist);
  else idx(e_dist > ref_electrode(2,1)) = []; 
  end

  m.electrode(ref+1).nodes = idx;
  m.electrode(ref+1).z_contact = 1; 
  m.gnd_node = m.electrode(ref+1).nodes(1);

  v_ = @(x) reshape(x,1,[]);
  
  % Update meas pattern
  meas_pattern = sparse(v_((1:nE)' * [1 1]), ... 
                           [1:nE (ref+1)*ones(1,nE)], ...
                            [ones(1,nE) -ones(1,nE)]);
  [m.stimulation.meas_pattern] = deal(meas_pattern); 

  for ii = 1:numel(m.stimulation) % resize sparse stim array
     m.stimulation(ii).stim_pattern(ref+1) = 0;
  end

  if chunk == 1 && any(named('-reference'))
    plot3(m.nodes(idx,1),m.nodes(idx,2),m.nodes(idx,3),'k.','MarkerSize',4)
    text(mean(m.nodes(idx,1)),mean(m.nodes(idx,2)),mean(m.nodes(idx,3)),'GROUND')
  end

end

return

%% use EIDORS.mk_image and fill in sigma values from the literature
function em = build_forward_model(fmdl,sigma)

% TODO - refactor this to load from the table in /source/ 

if ~exist('sigma','var')
    sigma(1).name = 'Interstitial';
    sigma(2).name = 'Prostate';      % 0.4243 S/m % http://dx.doi.org/10.1109/TBME.2007.897331
    sigma(3).name = 'PDMS';          % [materials properties] 
    sigma(4).name = 'Fascicle';      % 
    sigma(5).name = 'P_Fascicle';    % https://pubmed.ncbi.nlm.nih.gov/30507555/
    sigma(6).name = 'Epineurium';
    sigma(7).name = 'Blood';         % https://europepmc.org/article/pmc/pmc2914348

    sigma(1).value = 0.66;           % units of S/m
    sigma(2).value = 0.4243;
    sigma(3).value = 1e-12;
    sigma(4).value = [0.570, 0.088, 0.088]; 
    sigma(5).value = 8.7e-4;
    sigma(6).value = 0.083;          % [6.3] ? 
    sigma(7).value = 1.09; 

% Perineurium 0.0021
% Endoneurium
%  longitudinal 0.57
%  transverse   0.083
% Epineurium 0.083
% Saline 2.0
%
% from https://pubmed.ncbi.nlm.nih.gov/14278100/
% THE SPECIFIC IMPEDANCE OF THE DORSAL COLUMNS OF CAT
%
% from http://dx.doi.org/10.1016/0014-4886(65)90126-3
% Specific impedance of cerebral white matter. Nicholson (1965)
%     
%        Normal to fibers    Parallel to fibers
%        Resistive Reactive  Resistive Reactive
% 20 Hz  850 u 150  67 u 32   89 u 40  7 u 8
% 200 Hz 770 u 160  39 u 18   89 u 39  4 u 5
% 2 kHz  770 u 140  35 u 27   78 u 33  6 u 6
% 20 kHz 750 u 150  140 u 110 78 u 33  15 u 22
%
% Table 1: Conductivities used for modeling [12]: 
% AQ Choi, JK Cavanaugh, and DM Durand, "Selectivity of multiple-contact
% nerve cuff electrodes" IEEE Trans Biomed Eng (2001). 
    
end
if ~exist('fmdl','var'), em = sigma; return, end

em = mk_image(fmdl, 0.66);

em.elem_data(:,1,1:3,1:3) = 0;

for ii = 1:length(sigma)
    
  sel = strncmpi(em.fwd_model.object_name,sigma(ii).name, ...
                                 length(sigma(ii).name));
  sel = cat(1,em.fwd_model.object_id{sel});
  
  if ~any(sel), continue, end
  
  if numel(sigma(ii).value) == 3
    em.elem_data(sel,1,1,1) = sigma(ii).value(1);
    em.elem_data(sel,1,2,2) = sigma(ii).value(2);
    em.elem_data(sel,1,3,3) = sigma(ii).value(3);
  else
    em.elem_data(sel,1,1,1) = sigma(ii).value(1);
    em.elem_data(sel,1,2,2) = sigma(ii).value(1);
    em.elem_data(sel,1,3,3) = sigma(ii).value(1);
  end
end

% Check that every element has an assigned conductivity
done = any(em.elem_data(:,:,:),3);
if ~all(done)
    % warning('build_forward_model:unsetElements',
    error('%d of %d tets have undefined conductivities', ...
        sum(~done), length(done))
end

%% Quick mesh visualisation
function preview_eidors_mesh(m,do_animate)

clf, show_fem(m)
sel = cat(1,m.object_id{strncmp(m.object_name,'Elec',4)});
patch('Faces',m.elems(sel,:),'Vertices',m.nodes,'EdgeColor',[1 .4 .2], ...
                        'FaceColor','w','FaceAlpha',0.2,'EdgeAlpha',0.5)
h = get(gca,'Children');
h(end).EdgeAlpha = 0.2;
pause(0.1)

set(gca,'CameraPosition',get(gca,'CameraPosition') .* [-1 -1 1])
hold on, axis off

if nargin < 2, return, end
if ~do_animate, return, end

for ii = 3:numel(m.object_id)
    % be careful - there's two sets of indices for different objects which
    % aren't necessarially in the right order. 
    oid = m.object_id{ii};
    h(1).Faces = m.elems(oid,:); 
    title(m.object_name{ii})
    ginput(1);
end

clear oid i h b t

%% Convert v_extracellular to structure as output by "coalesce"
function EM = convert_stim2fascicles(EM)
% see https://en.wikipedia.org/wiki/Reciprocity_(electromagnetism)#Reciprocity_for_electrical_networks               

  nF = sum(contains(EM.model.object_name,'Fasc')) - ...
       sum(contains(EM.model.object_name,'P_Fasc'));
    
  fasc_ = @(n) sprintf('Fascicle%d',n);
  
  % EM.v_extracellular is a list, (nodes x electrodes) of the sequential
  % bipolar stimulus-induced extracellular potential, for ALL nodes 
  % (not just fascicle nodes). 

  % a sensitivity file is broken up by fascicle, which each fascicle
  % containing a NODE index and a list of per-node values 
  
  % EM.model.object_id is a list of ELEMENT ids which make up each object
  % in the EIDORS model. Gather each NODE in that set of elements, then use
  % that to emulate the EM.FascicleN.pot/idx 
    
  for ff = 1:nF
        
    sel = strcmpi(EM.model.object_name,fasc_(ff));
    if sum(sel) == 0 && FLAG_multimesh 
      sel = strcmpi(EM.model.object_name,fasc_(1));      
    end
    if sum(sel) ~= 1, error('Object %s not found in EM.model.object_name',fasc_(ff)); end
    idx = unique(EM.model.elems(EM.model.object_id{sel},:));
    
    EM.(fasc_(ff)).idx = idx'; 
    EM.(fasc_(ff)).pot = EM.v_extracellular(idx,:)';
  end  
  
  EM.utils.x_ = '@(i) out.model.nodes(i,1)';
  EM.utils.y_ = '@(i) out.model.nodes(i,2)';
  EM.utils.z_ = '@(i) out.model.nodes(i,3)';
  EM.utils.fac_ = '@(i) out.(sprintf(''Fascicle%d'',i)).idx'; 

return
