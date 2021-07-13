
function axons = axon_sfap(varargin)
% models.axon_sfap computes single-fibre action potential magnitudes for 
%   every axon in an input population (generated using axon_population). 
%   models.axon_sfap requires a sensitivity field (generated using 
%   nerve_anatomy and membrane current profiles (generated using 
%   membrane_currents], and is structurally similar to 
%   models.nerve_recording. 
% 
% In addition to the '-name', value syntax, models.axon_sfap accepts an 
%   options structure; the fields of `options.sfap` are treated as 
%   additional input arguments (see tools.opts_to_args). 
% 
% The output SFAP wave is in units of:
% SFAP (uV) =  I (nA) / [ 4 pi sigma (ohm/m) distance (mm)]
% The mm and the nA cancel to yield �V, so no further conversion is needed.
% 
% Core Syntax: 
%  -file 'sensitivity*.mat' : specify spatial sensitivity function
%  -axons 'axons*.mat' : specify axon population for SFAP calculation, 
%                        defaults to newest axons file in `axons~/`
%  -currents [folder] : specify root folder for membrane currents.
%  -out : specify output folder location (default: `sfap~/`)
% 
% Additional options: 
% -axon-t : set 3D axon trajectories (not yet implemented, may interfere 
%           with -axons?)
% -no-parallel : run in debug mode (without using the parallel pool)
% -recenter-peak [new_peak_x]: shift simulated sensitivity peaks in space. 
%                              Useful for the situation in which many 
%                              electrodes of different designs were 
%                              simulated in a single mesh. 
%                              If numel(new_peak_x) > 1, modulo. 
% -xy : replace axon positions with a homogenous grid of axons. This is 
%       actually particularly useful with models.axon_sfap. The grid is 
%       based on the median axon for each class; this can be overridden 
%       with the additional arguments:  -ad [myelinated axon diameter],
%                                       -gr [g-ratio], 
%                                       -cd [unmyelinated axon diameter]. 
%       If the argument -xy-eidors is passed, the axon grid ignores the 
%       nerve profile in the axons file and instead uses a grid based on 
%       the fascicle geometry in the sensitivity file. -xy-quick generates 
%       a smaller 1D sample and -xy-resol [n_axons] sets grid resolution. 
% -distortion : apply a systemic distortion to the recorded velocities, 
%               temporal current width, and/or spatial current width 
%               (not tested in a long time). -debug-distortion illustrates.
% -ref [value] : reference x-axis value for raster spike times
% -fs [rate, default 30 kHz]   : set sampling rate for recorded SFAP
% -time [value, default 20 ms] : set time ROI window for SFAP
% -q,-quiet : suppress most output to console. 
% -debug-units : set up the calculation and display axon trajectories 
%                against the fascicle, useful for debugging units issues. 
%                Synonym for -preview
% 
% Last updated 23-June-2021 CDE v0.3


varargin = tools.opts_to_args(varargin,'sfap');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('-q')), printf = @(varargin) []; else printf = @fprintf; end

printf('Running models.%s ... \n', mfilename);
if isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP>
    save_default_options ('-mat-binary'), 
end

[EM,AX] = tools.parse_arguments(varargin, 'LOAD', ...
                              'eidors','eidors~/sens*.mat', 'axons');
                          
output_file = EM.filename; 
if any(named('-out')), output_file = get_('-out'); end
if isa(output_file,'function_handle'), output_file = output_file(); end
if ~any(ismember(output_file,'/\')), % if not already a path ... 
  output_file = tools.file(['sfap~/' output_file '.mat'],'-make');
end

while exist(output_file,'file')
  output_file = [output_file(1:end-4) '_NEW.mat'];
end

if any(named('-axon-t')), EM.info.AxonTrajectory = get_('-axon-t'); end

printf('Results will be saved to %s \n', tools.file('TT',output_file));

opts = struct;
opts.axons_folder = AX.folder; 
opts.verbose = ~any(named('-q'));
opts.use_parallel = ~any(named('-no-p'));

pop = AX.pop; 
nerve = AX.nerve; %#ok<NASGU>

%% Get sensitivity functions 

for f = fieldnames(EM.utils)' % Create utility local functions
  EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');
  eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
end

fascicles = fieldnames(EM); % "Fascicle1" not guaranteed to exist
fascicles(~strncmp(fascicles,'Fascicle',8)) = [];

nE = size(EM.(fascicles{1}).pot,1);
nF = sum(strncmp(fieldnames(EM),'Fascicle',6)) ; % size(nerve.fascicles,3); 

sensitivity = cell(nE,nF); % i2v (uV/uA) as f(x,y,z) 

for ff = 1:nF % get sensitivity function 
  sensitivity{1,ff} = scatteredInterpolant( ... 
                       z_(fac_(ff)), y_(fac_(ff)), x_(fac_(ff)), ...
                       EM.(sprintf('Fascicle%d',ff)).pot(1,:)', ...
                       'natural','none'); % units of mm ? 
  for ee = 2:nE % copy to speed up 
     sensitivity{ee,ff} = sensitivity{1,ff};
     sensitivity{ee,ff}.Values = EM.(sprintf('Fascicle%d',ff)).pot(ee,:)';
 end
end

% Translate sensitivity peaks if specified
if any(named('-recenter-peak'))
  if verbose, disp('Translating sensitivity peaks...'), end
  sensitivity = translate_I2V_peaks(sensitivity, get_('-recenter-peak'));
end

% If -grid-XY replace axon positions with XY grid
if any(named('-xy')),
  %%
  if verbose, disp('using XY grid'), end
  pop = make_axon_xy_grid(EM,pop,AX.nerve); 
  nF = max(pop.unmyelinated_fascicle);
end

clear ee ff x_ y_ z_ fac_ f sel ok

%% Get Distortion if specified

if any(named('-dis')),
  opts.distortion = get_('-dis'); 
  if size(opts.distortion,1) == numel(pop)
    opts.distorion_by_type = opts.distortion; 
  end
  % options.distortion. [ v z t ] is read as a matrix or structure by
  % models.spike_to_wave.
  
  if any(named('-debug-dis')), opts.plot_image = true; end  
end

% default position for spiketime: a nominal spike-time of "0" 
%   means the AP is at z=0 at time 0 (axons oriented along
%   z-axis and propegate in the +z direction) 
z_ref = 0; 
if any(named('-ref')) z_ref = get_('-ref'); end
if ischar(z_ref), z_ref = str2double(regexp(z_ref,'\d+','match','once'));
    z_ref = median(EM.model.nodes(EM.model.electrode(z_ref).nodes,1)); 
end

fs = 30; % 24.414; % kHz 
if any(named('-sample-rate')), fs = get_('-sample-rate'); 
elseif any(named('-fs')),      fs = get_('-fs'); 
end

time_span = 20;
if any(named('-time')),  time_span = get_('-time'); end

%% Compute sensitivity at each fibre xy and take PCA

sam = struct; 

for myelin = 0:1 % for each axon class 
    
  % I think this falls over once different myelinated axons have
  % different ultrastructural relationships. 
  
  sel = find([pop.myelinated] == myelin); 
  [subtype,group_id,axon_xy] = get_axon_xy(EM, pop(sel));
  
  if isempty(group_id), continue, end  
  printf('Computing sensitivity PCA (myelin=%d) ... \n',myelin);
    
  nG = max(group_id); 
  nF = max(cat(1,pop(sel).fascicle));

  % axon_xy is a {nG x nF} cell array of axon positons. By default, 
  % each element is a [nCells x 2] list of intercepts in the xy plane.
  % If EM.info.AxonTrajectory or EM.info.FascicleTrajectory is defined,
  % each element is instead a [nCells x 3 x nPoints] trajectory for
  % each cell. 
  
  if any(named('-debug-u')) || any(named('-preview'))
    %% Debug this, is a consistent problem ... 
    tools.setupEIDORS ;
    clf, plots.view_mesh(EM.model,'-fasc');
    hold on, h = get(gca,'Children');
    set(h,'EdgeAlpha',0.1)

    for ii = 1:nG*nF
      for uu = 1:size(axon_xy{ii},1)
        xyz = permute(axon_xy{ii}(uu,:,:), [2 3 1]);
        if size(axon_xy{ii},2) == 2
          plot3(0*xyz(1,:), xyz(2,:),xyz(1,:),'o','MarkerSize',4, ...
              'Color',[1 0 0 0.2],'Clipping','off') 
        else
          plot3(xyz(1,:), xyz(2,:),xyz(3,:), ...
              'Color',[1 0 0 0.2],'Clipping','off') 
        end
        hold on
      end
    end

    if size(axon_xy{ii},2) == 2
      view([1 0 0]), axis([min(xlim) max(xlim) xyz(2)-[0.2 -0.2] ...
                                               xyz(1)-[0.2 -0.2]])
    end
    %%
    return
  end
  
  i2v_pca = cell(size(axon_xy)); 
    
  for ii = 1:(nG*nF) % The parfor moved to inside get_i2v_pca  
    
    [gg,ff] = ind2sub([nG nF],ii);
    i2v_pca{ii} = get_i2v_pca(sensitivity(:,ff),axon_xy{ii}, ... 
                          subtype(gg),pop(sel(1)).axon_model,opts);  %%#ok<PFBNS>
    assert(isempty(i2v_pca{ii}.weight) || any(i2v_pca{ii}.profile(:)), ...
           'Computed 0 or NaN PCA profile')
  end
  
  i2v_pca = reshape([i2v_pca{:}],nG,nF);  
  the = @(v,x) arrayfun(@(g) median(v(g==x)), (1:nG)');

  pop_id = arrayfun(@(x) x*ones(size(pop(x).size_sample)), sel,'unif',0);
  
  switch myelin
    case 0 % all unmyelinated fibres 
      sam.C_populations = cat(1,pop_id{:}); 
      sam.C_index = group_id;
      sam.C_fascicle = cat(1,pop(sel).fascicle);
      sam.C_groups = subtype;      
      sam.C_diam = the(cat(1,pop(sel).fibre_diam), group_id); 
      sam.C_PCA_sensitivity = i2v_pca;
    case 1 % all myelinated fibres
      sam.A_populations = cat(1,pop_id{:}); 
      sam.A_index = group_id;
      sam.A_fascicle = cat(1,pop(sel).fascicle);
      sam.A_groups = subtype;      
      sam.A_diam = the(cat(1,pop(sel).fibre_diam), group_id); 
      sam.A_g_ratio = the(cat(1,pop(sel).g_ratio), group_id); 
      sam.A_PCA_sensitivity = i2v_pca;
  end  
end

clear myelin sel group_id subtype ff gg the i2v_pca
% "sam" is now set up and equiped with A_ and C_PCA_sensitivity

%% Compute intracellular-to-extracellular relationship

axons = []; 
warn_once = true; 

for ty = 1:numel(pop)
  time = 0:(1/fs):(time_span);
  time = [-fliplr(time) time(2:end)]; % ms

  opts.filename = '';
  opts.reference_z = z_ref; 
  opts.class = pop(ty).axon_model;
  opts.padding = 2;
  opts.efferent = ~pop(ty).afferent;
  
  [subtype,index,axon_xy] = get_axon_xy(EM, pop(ty));
  axon_xy(all(cellfun(@isempty,axon_xy),2),:) = []; 
  
  if pop(ty).myelinated      
      subtype = find(ismember(sam.A_groups,subtype));      
      opts.sensitivity_pca = sam.A_PCA_sensitivity(subtype,:); 
  else
      subtype = find(ismember(sam.C_groups,subtype));      
      opts.sensitivity_pca = sam.C_PCA_sensitivity(subtype,:);
  end
  
  nG = numel(subtype);
  V = cell(nG,nF);
  
  cache_path = tools.cache('path');

  %% CORE parfor loop 
  
  if nG*nF == 1 || ~opts.use_parallel
    for ii = 1:(nG*nF)
      V{ii} = parfun_unpack(cache_path, ...
                     @(g,f) pca_to_wave(subtype(g),time, ...
                                        opts.sensitivity_pca(g,f),...
                                        axon_xy{g,f},opts), ...
                           [nG nF],ii);
    end
  elseif isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP> % octave parallel
    if isempty(which('pararrayfun')), pkg load parallel; end 
    V = pararrayfun(nproc-1, @(a) parfun_unpack(cache_path, @(g,f) ...
                                   pca_to_wave(subtype(g),time, ...
                                        opts.sensitivity_pca(g,f),...
                                        axon_xy{g,f},opts), [nG nF], a), ...
                    1:(nG*nF), 'UniformOutput', false);
  else    
    parfor ii = 1:(nG*nF)  
      V{ii} = parfun_unpack(cache_path, ...
                     @(g,f) pca_to_wave(subtype(g),time, ...
                                        opts.sensitivity_pca(g,f),...
                                        axon_xy{g,f},opts), ...
                           [nG nF],ii); %#ok<PFBNS>
    end
  end
  
  %% Format output to something useful
  
  this = struct;
  
  this.class = pop(ty).axon_model; 
  if opts.efferent, this.class = [ this.class ' (efferent)'];
  else              this.class = [ this.class ' (afferent)'];
  end
  
  this.time = time;
  this.axon.xy = pop(ty).axon_xy; 
  this.axon.fascicle = pop(ty).fascicle;
  this.axon.subtype = index;
  this.axon.subtype_index = subtype;
  this.sensitivity = opts.sensitivity_pca; 
  
  for ii = 1:nG*nF
      [gg,ff] = ind2sub([nG nF],ii); 
      if isempty(this.sensitivity(ii).weight), continue, end      
      if pop(ty).myelinated
           sel = sam.A_populations == ty; 
           sel = sel(sam.A_index == subtype(gg) & sam.A_fascicle == ff);
      else sel = sam.C_populations == ty; 
           sel = sel(sam.C_index == subtype(gg) & sam.C_fascicle == ff);
      end      
      this.sensitivity(ii).info = 'computed using PCA()';     
      this.sensitivity(ii).weight(~sel,:) = []; 
      this.sensitivity(ii).xy = axon_xy{gg,ff}; 
  end
        
  if pop(ty).myelinated
       this.axon.diameter = sam.A_diam(subtype);
       this.axon.g_ratio  = sam.A_g_ratio(subtype);
  else this.axon.diameter = sam.C_diam(subtype);
  end
  
  if numel(V) == 1, this.component_SFAP = V{1};
  else this.component_SFAP = V;
  end
  
  this.units.time = 'ms';
  this.units.xy = 'mm';
  this.units.axon_subtype = '1 .. 24 from models.axon_population';
  this.units.sensitivity_profile = 'length x electrode x nK, uV/uA';
  this.units.sensitivity_weight  = 'nAxon x nK, unitless';
  this.units.sensitivity_latent  = '% variance explained';
  this.units.component_SFAP = 'time x electrode x nK, uV';
  
  % This check prevents from having to open another instance  
  su_ = @(x) reshape(x,[],1);
  su_ = @(x) [sqrt(nanmean(su_(x).^2)) nanmin(su_(x)) nanmax(su_(x)) ];  
  summary = su_(cat(3,V{:}));
  printf(' %s: RMS %g, range %g-%g\n', this.class, summary ); 
  if warn_once && all(summary == 0)
      warning('ViNERS:possibleUnitsFail','WAVE was all-zeros, try calling with -unit-um or -debug-units')
      warn_once = false; 
  end
  
  if numel(V) == 1, 
    nK = numel(this.sensitivity.latent);
    this.axon_SFAP = reshape(V{1},[],nK) * this.sensitivity.weight' ;
    this.axon_SFAP = reshape(this.axon_SFAP, ...
                             size(this.component_SFAP,1), ...
                             size(this.component_SFAP,2), []);
    this.units.axon_SFAP = 'time x electrode x nAxon, uV';      
  else
    
    this.axon_SFAP = zeros(size(V{1},1),nE,numel(this.axon.subtype));
    this.units.axon_SFAP = 'time x electrode x nAxon, uV';
    
    for ii = 1:(nG*nF)
      [gg,ff] = ind2sub([nG nF],ii); 
      nK = numel(this.sensitivity(ii).latent);

      grp_sfap = reshape(V{ii},[],nK) * this.sensitivity(ii).weight' ;
      grp_sfap = reshape(grp_sfap, size(V{ii},1), size(V{ii},2), []);
      this.axon_SFAP(:,:,index==gg & pop(ty).fascicle==ff ) = grp_sfap;
    end
  end
  
  if isempty(axons), axons = this;
  else axons(ty) = this; %#ok<AGROW>
  end

end

%% Save output to sub-xx/sfap/(eidors_filename).mat

clear axon_xy xy grp fid index opts nG z_ref V ty ff 
clear mdl myelin ii is_hash

eidors_file = EM.filename;
nerve = AX.nerve;

printf('Saving %s\n ', tools.file('TT',output_file))
save(output_file,'axons','eidors_file','nerve','pop','sam');

if nargout == 0, clear, end
return

%% Move sensitivity peak (for electrode array property explorations) 
function sensitivity = translate_I2V_peaks(sensitivity, z_REF)
%%

if numel(z_REF) < nE, 
  z_REF = repmat(z_REF(:), [ceil(nE/numel(z_REF)) 1]); 
end

color = 1-summer(nE);
figure(3), clf

for ff = 1:nF % for each fascicle 
  xy = [mean(z_(fac_(ff))) mean(y_(fac_(ff)))];
  z = linspace(min(x_(fac_(ff))), max(x_(fac_(ff))),1001);
  subplot(nF,1,ff), hold on
  for ee = 1:nE

    pk0 = sensitivity{ee,ff}([xy 0] + z'*[0 0 1]); 
    [~,id] = nanmax(pk0);   
    sensitivity{ee,ff}.Points(:,3) = sensitivity{ee,ff}.Points(:,3) ...
                                     - z(id) + z_REF(ee);
    pk1 = sensitivity{ee,ff}([xy 0] + z'*[0 0 1]);    
    plot(z,pk0,':','Color',color(ee,:),'LineWidth',1.2)
    plot(z,pk1,'-','Color',color(ee,:),'LineWidth',1.2)
    % pause(0.02)
  end
  tools.tidyPlot, pause(0.05)
end 

return

%% Throw out xy coordinates from axons~/axons.mat and use an XY grid
function pop = make_axon_xy_grid(EM,pop,nerve)

named = evalin('caller','named');
get_ = evalin('caller','get_');

if ~isfield(nerve,'fascicles'), nerve.fascicles = nerve.outline; end 

use_F_anatomy = ~any(named('-xy-e')); % XY-em or xy-eidors
if use_F_anatomy, nF = size(nerve.fascicles,3);
else nF = sum(strncmpi(EM.model.object_name,'Fascicle',6));
end

resol = 11; 
use_1d_sample = any(named('-xy-q')); % XY-grid or XY-quick
if any(named('-xy-c')), resol = get_('-xy-c'); end  
if any(named('-xy-r')), resol = get_('-xy-r'); end  

%%
if ~use_F_anatomy
  fac_ = @(n) unique( EM.model.elems( ...
                      EM.model.object_id{strcmp(EM.model.object_name, ...
                                        sprintf('Fascicle%d',n))}, : ));  
  % x_ = @(n) EM.model.nodes((n),1); 
  y_ = @(n) EM.model.nodes((n),2); 
  z_ = @(n) EM.model.nodes((n),3); 
end

% source = pop;

for ty = 1:numel(pop)
            
  pop(ty).fibre_diam = median(pop(ty).fibre_diam);
  if pop(ty).myelinated

    if any(named('-ad')), pop(ty).fibre_diam = get_('-ad'); end
    if any(named('-ag')), pop(ty).g_ratio = get_('-ag'); end
    if any(named('-gr')), pop(ty).g_ratio = get_('-gr'); end
      
    pop(ty).g_ratio = median(pop(ty).g_ratio);
    pop(ty).axon_diam = pop(ty).fibre_diam .* pop(ty).g_ratio;
    
    if ~any(named('-q')), 
      fprintf('[%02d | %s] %s, d=%0.4f, g=%0.4f\n', ty, ...
            pop(ty).axon_model,'Simulating grid of axons', ...
            pop(ty).fibre_diam, pop(ty).g_ratio); 
    end        
  elseif ~any(named('-q')), 
      if any(named('-cd')), pop(ty).fibre_diam = get_('-cd'); end
      fprintf('[%02d | %s] %s, d=%0.4f\n', ty, ...
            pop(ty).axon_model,'Simulating grid of axons', ...
            pop(ty).fibre_diam); 
  end
  pop(ty).size_sample = 0; 
  
  if ty > 1      
    pop(ty).axon_xy = pop(1).axon_xy; 
    pop(ty).fascicle = pop(1).fascicle; 
    continue
  end

  %% Get XY positions for each fascicle 
  pop(ty).axon_xy = []; 
  pop(ty).fascicle = []; 
  
  for ff = 1:nF
  
    if use_F_anatomy, 
      xy_fac = unique(nerve.fascicles(:,:,ff),'rows','stable'); 
    else
      xy_fac = [z_(fac_(ff)) y_(fac_(ff))];
      sel = convhull(xy_fac(:,1),xy_fac(:,2));
      if numel(sel) > 100, sel = sel(round(linspace(1,end,101))); end
      xy_fac = xy_fac(sel([1:end 1]),:);
    end
  
    if numel(resol) > 1 % explicit axon coordinates 

      gx = resol(:,1); 
      gy = resol(:,2); 
      ok = inpolygon(gx(:),gy(:),xy_fac(:,1),xy_fac(:,2));
      
    elseif use_1d_sample % -xy-quick : line sample from min-y to max-y
    
      [~,idx] = min(xy_fac(:,2));
      [~,idx(2)] = max(xy_fac(:,2));
    
      gx = linspace(xy_fac(idx(1),1),xy_fac(idx(2),1),resol+1);
      gy = linspace(xy_fac(idx(1),1),xy_fac(idx(2),1),resol+1);
      gx = conv(gx,[1 1]/2,'valid')'; 
      gy = conv(gy,[1 1]/2,'valid')'; 
      ok = true(size(gx));
    
    else % -xy-grid
    
      xy0 = [median(xy_fac) min(xy_fac) max(xy_fac)];
      xy0(3:4) = max(abs(xy0(1:2)-xy0([3 4; 5 6])),[],1);
      [gx,gy] = meshgrid(xy0(1)-xy0(3)*linspace(-1,1,resol), ...
                         xy0(2)-xy0(4)*linspace(-1,1,resol));
      ok = inpolygon(gx(:),gy(:),xy_fac(:,1),xy_fac(:,2));
    end
  
    pop(ty).axon_xy = [pop(ty).axon_xy; gx(ok) gy(ok)];
    pop(ty).fascicle = [pop(ty).fascicle; 0*gx(ok) + ff];
  end
end
  
n_axons = size(pop(1).fascicle); 
for ty = 1:numel(pop)
  for f = {'fibre_diam','axon_diam','g_ratio'}
    if isempty(pop(ty).(f{1})), continue, end
    pop(ty).(f{1}) = pop(ty).(f{1})(1) * ones(n_axons);
  end
end

return

%% Loop over axons to get the sensitivity curve (contains a call to parfor)
function [p,i2v] = get_i2v_pca(sensor,xy,gg,name,opts)

if isfield(opts,'axons_folder')
     d = sprintf('%s%s/n%03d_out.mat',opts.axons_folder,name,gg);
else d = tools.file(sprintf('axons~/%s/n%03d_out.mat',name,gg));
end

d = load(d,'index','length'); 
xy_list = xy; 

%% WAVE = �V = (nA) / (4*pi*sigma) (ohm.m) / distance (mm) 
% the mm and the nA cancel to yield �V, no conversion necessary

[Z,~] = sort(d.length);
nZ = numel(Z);
nE = numel(sensor); 
nA = size(xy_list,1);
i2v = zeros(nZ,nE,nA); 

if ndims(xy_list) == 3 % 3D axon trajectories, this whole approach might fall over...
  error TODO_check_this_code
  xy = permute(xy,[3 2 1]); %#ok<UNRCH>
  u = sqrt(sum(diff(xy).^2,2)); 
  u = cumsum(conv(u([1 1:end end]),[1;1]/2,'valid'));
  [~,z0] = min(abs(xy(:,1))); u = u-u(z0);

  Zq = fliplr(interp1(u,xy,Z));
  i2v = sensor{ee}(Zq(:,1),Zq(:,2),Zq(:,3));

else
  if ~opts.use_parallel % any(named('-no-p'))
    for ii = 1:size(xy_list,1)                    
      i2v(:,:,ii) = get_i2v(sensor,xy_list(ii,:),Z); 
    end   
  elseif isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP> % octave parallel
    if isempty(which('pararrayfun')), pkg load parallel; end      
    i2v = pararrayfun(nproc-1, @(a) get_i2v(sensor,xy_list(a,:),Z), ...
                                    1:size(xy_list,1),'UniformOutput',0);
    i2v = cat(3,i2v{:}); % convert to struct array      
    error TODO_validate_this_octave_code
  else 
    parfor ii = 1:size(xy_list,1)
      i2v(:,:,ii) = get_i2v(sensor,xy_list(ii,:),Z); 
    end
  end
end

%%
if isempty(i2v) % nothing to simulate? 
  
 p.weight = zeros(0,1);
 p.latent = 100;
 p.profile = zeros(nZ,nE,1);
 return

elseif all(i2v(:) == 0)
  error(['The axon population did not intersect the fascicle geometry, ' ... 
         'please check your units and try again or call using -xy-eidors'])
end

i2v = permute( i2v, [3 1 2] ); % axonID, elec, z       
[p.profile,p.weight,p.latent] = pca(i2v(:,:),'centered',false);

scale = max(p.weight,[],1); 

null = (scale == 0); % an axon outside of fascicle?
scale(null) = []; 
p.profile(:,null) = []; 
p.weight(:,null) = []; 
p.latent(null) = []; 

p.profile = p.profile .* scale;
p.weight  = p.weight  ./ scale;

assert(~any(isnan(p.weight(:))),'NaN weight computed, please check')

nK = find(p.latent / sum(p.latent) < 0.005,1);

if isempty(nK), nK = 1; end

p.weight = p.weight(:,1:nK);
p.latent = 100 * p.latent(1:nK) / sum(p.latent);
p.profile = reshape(p.profile(:,1:nK),[],nE,nK);


return

function i2v = get_i2v(sensor,xy_list,Z)

  nZ = numel(Z);
  nE = numel(sensor); 
  i2v = zeros(nZ,nE); 
  
  xy = xy_list;
  for ee = 1:nE
      
    pot = sensor{ee}(0*Z+xy(1), 0*Z+xy(2), Z); %%#ok<PFBNS>
    if all(isnan(pot(Z >= 0))) && ~all(isnan(pot(Z <= 0)))    
      % To save time and space, only the right half (+z) was simulated.
      %   Use symmetry to get the other half of the sensitivity curve. 
      pot(Z<0) = sensor{end-ee+1}(0*Z(Z<0)+xy(1), 0*Z(Z<0)+xy(2), -Z(Z<0));
    end

    ok = ~isnan(pot); % Patch over discontinuities
    if mean(ok) < 0.5, i2v(:,ee) = 0; continue, end
    if ~all(ok),
      pot(~ok) = interp1(Z(ok),pot(ok),Z(~ok),'linear','extrap');
    end

    i2v(:,ee) = pot;
  end
  
% WAVE = uV = (nA) / (4*pi*sigma) (ohm.m) / distance (mm)
% the mm and the nA cancel to yield uV, no conversion necessary


%% Get a cell array of axon XY for the contents of the population
function [sample_ids,gid,xy] = get_axon_xy(EM,pop)

[sample_ids,~,gid] = unique(cat(1,pop.size_sample)); 
% 
%   % axon_xy is a {nG x nF} cell array of axon positons. By default, 
%   % each element is a [nCells x 2] list of intercepts in the xy plane.
%   % If EM.info.AxonTrajectory or EM.info.FascicleTrajectory is defined,
%   % each element is instead a [nCells x 3 x nPoints] trajectory for
%   % each cell. 
%   

nF = max(cat(1,pop.fascicle)); 
nG = numel(sample_ids);

fid = cat(1,pop.fascicle); 
source_xy = cat(1,pop.axon_xy);

xy = cell(nG,nF); 
nerve = evalin('caller','nerve'); 
      
if tools.from_trajectory(EM,nerve)
  for ff = 1:nF % For each fascicle
    xy(:,ff) = tools.from_trajectory(EM,nerve,source_xy(fid == ff,:), ...
                                            '-f',ff,'-g',gid );
  end
else
  for ii = 1:(nG*nF) 
    [gg,ff] = ind2sub([nG nF],ii); 
    xy{ii} = source_xy(gid == gg & fid == ff,:); 
  end % setup g_xy    
end
  
%% Compute PCA component contribution for any given spike 
function wave = parfun_unpack(cache_path,fun,nGnF,ii)

  tools.cache('set',cache_path)
  [gg,ff] = ind2sub(nGnF,ii);
  wave = fun(gg,ff); 
%  V{ii} = pca_to_wave(index(gg),time,opts.sensitivity_pca(gg,ff),...
%                                   axon_xy{ii},opts); %%#ok<PFBNS>


function wave = pca_to_wave(g_id,time,sensors,~,options)
% based on models.spike_to_wave which retreives the saved membrane current 
%  for a specified axon (g_id) and simulates the extracellular potential. 
% 
% g_id = axon template index
% time = simulation time-points
% sensors = sensitivity (pca structure)
% g_xy = xy co-ordinates of axons (units to match sensors) 
% % if g_xy is 2D, assume straight fascicles; 

% options
%   .raster = models.random_raster % what spike pattern to simulate?
% 
%   .axons_folder / .class % what membrane current profiles to use? 
% 
%   .reference_z = 0 % what point in space are spikes aligned to?
%   .padding = 1     % left/right border padding on sensitivity function
%   .evoked_potential = [] % ? 
%   .efferent = 0 % flip propegation direction
%   .distortion = [1 1 1] % edit CV/Z-width/T-width
%   .get_all = 0 % return partial sum of the wave for each spike in _spk1
%   .plot_image = 0 % Make graphical output
% 
% TODO - make this code faster. Almost half of all model time is spent in
%          this function. 
%
% version 0.5 - 10 July 2020 Calvin Eiber
% version 0.4 - 27 May 2020  Calvin Eiber) 

if ~exist('options','var'), options = struct; end
if ~isfield(options,'reference_z'), options.reference_z = 0; end
if ~isfield(options,'padding'), options.padding = 1; end
if ~isfield(options,'evoked_potential'), options.evoked_potential = []; end
if ~isfield(options,'plot_image'), options.plot_image = false; end

nE = size(sensors(1).profile,2);
nK = size(sensors(1).profile,3); 

wave = zeros(length(time),nE);
spike_time = 0;

if isfield(options,'axons_folder')
     d = sprintf('%s%s/n%03d_out.mat',options.axons_folder,options.class,g_id);
elseif isfield(options,'class')
     d = tools.file(sprintf('axons~/%s/n%03d_out.mat',options.class,g_id));
else d = sprintf('%sn%03d_out.mat',tools.cache('PATH'),g_id);
end

d = load(d); 

if isfield(options,'distortion')
  distortion_factor = options.distortion; 
  if isstruct(distortion_factor),     
      distortion_factor = [options.distortion.v ...
                           options.distortion.z ...
                           options.distortion.t];
  end
else distortion_factor = [1 1 1];
end
if isfield(options,'efferent') && options.efferent
   distortion_factor = [-1 1 -1] .* distortion_factor;
end
    
nZ = (numel(d.index)-1);
nP = (numel(d.length)-1) / nZ + 1;    

[~,t0] = max(d.V_example(:,2)); 
[~,p0] = min(abs(d.length(1:nP) - options.reference_z));

time = reshape(time,[],1);
% d.dt_dx is the delta-ms per length sample 

[Z,~] = sort(d.length);

assert(numel(Z) == size(sensors(1).profile,1),... 
      'mismatch between PCA I2V (nZ=%d) and n%03d_out.mat (nZ=%d)', ... 
         size(sensors(1).profile,1), g_id, numel(Z))

dx = mean(diff(time)) / d.dt_dx; % integration scale for cumsum() 


for kk = 1:nK % for each PCA component ... 
  % printInfo('Group %d [u%03d]', gg,ii)
  
  c = d.I_membrane(:,2:end) ; % current in nA -> uA
  % c(c<0) = c(c<0) / sum(c(c<0)) * -sum(c(c>0)); % Balance currents
  v = d.dt_dx;
  
%   if isempty(distortion_factor), t = (spike_time(ii)-d.time+d.time(t0));
%     current = @(n) interp1(t, c, time-n*v,'Linear',0); % at node n, current at each time
%   else t = spike_time(ii)-(d.time-d.time(t0))/distortion_factor(3);
%     current = @(n) interp1(t,c,time-(n*v*distortion_factor(1)),'Linear',0);  
%   end

  if isempty(distortion_factor), t = spike_time + (d.time-d.time(t0));
    current = @(n) interp1(t, c, time-n*v,'Linear',0); % at node n, current at each time
  else t = spike_time + (d.time-d.time(t0))/distortion_factor(3);
    current = @(n) interp1(t,c,time-(n*v/distortion_factor(1)),'Linear',0);  
  end
  if ~isempty(options.evoked_potential), 
    poststimulus = 1./(1+exp((spike_time-time+options.evoked_potential) ...
                                / options.raster_opts.tau1 * 3));
    current = @(n) current(n) .* poststimulus;
  end
  
  clear c v t

  %%  
  if options.plot_image
    %% DEBUG image to check prop velocity and spike timing
    % this is basically a panel of the figure generated by
    % distortion_visualisation
    distortion_image = [];
    ee = 2; 
    
    i2v = sensors(1).profile(:,ee,kk); 
    
    % tx = tic;
    for pp = 1:nP-1
      i2v_seg = i2v( (pp-1)*nZ+(1:nZ) );
      distortion_image = [distortion_image current(pp-p0) * i2v_seg];
    end
    % tx = toc(tx);
  
    clf, % subplot(3,1,[1 2])
    imagesc(time, Z, log10(abs(distortion_image))')
    hold on, plot(0,options.reference_z,'r+','markersize',30)
    xlabel('time, ms'), ylabel('Z, mm')    
    
    plot(time,sum(distortion_image,2),'-w','LineWidth',1.1)
    % error('Generated example distortion image')
    return
  end
  
  %%
  
  i2v = sensors.profile(:,:,kk); % z-coord, elec
  if mean(i2v(:) == 0) > 0.5, continue, end    
    
  %%
  V = 0*wave(:,:,1); 
  balance_ = @(x) -x(x<0).*sum(x(x>0))./sum(x(x<0)); 
    
  if ~isempty(options.padding)

    C0 = current(0);  C0(C0<0) = balance_(C0);
    i2v_seg = ones(nZ,1)*i2v(1,:);
    
    edge_effect = cumsum(C0 * i2v_seg,'reverse');
    adi = mean(abs(diff(edge_effect([1 1:end],:))),2);     
    aoi = (adi > 1.5 * median(adi(adi>0)));
    nPad = ceil(sum(aoi) * dx * options.padding);
    vPad = 0*V;

    for pp = -nPad:0, % LHS padding
      C = current(pp-p0);
      % C(C<0) = balance_(C); 
      vPad = vPad + (C * i2v_seg);
    end
    
    i2v_seg = ones(nZ,1) * i2v(end,:);
    for pp = nP+(1:nPad+1), % RHS padding
      C = current(pp-p0);
      % C(C<0) = balance_(C); 
      vPad = vPad + (C * i2v_seg);
    end
    
    V = V + vPad;

  else nPad = 0; vPad = 0; 
  end

  
  %% Spatial summation of fields as the spike passes over the electrode
  
  for pp = 1:nP-1 % insert points
      C = current(pp-p0);
      C(C<0) = balance_(C); 
      i2v_seg = i2v((pp-1)*nZ+(1:nZ),:);
      V = V + C * i2v_seg;
  end

  % last segment - extrapolate i2v
  i2v_seg = i2v((pp-1)*nZ+(1:nZ),:) - i2v(pp*nZ-nZ+1,:) + i2v(pp*nZ+1,:);
  V = V + current(pp-p0+1) * i2v_seg;

  %% Apply some corrections to get rid of edge artifacts 
  % (not sure how effective this lot actually is)
  
  % C0 = current(0-nPad-p0-1.5);  % C0(C0<0) = balance_(C0);
  % C1 = current(nP+nPad+2.5-p0); % C1(C1<0) = balance_(C1);
  % edge_effect = cumsum(C1 * V2I_seg) + cumsum(C0 * i2v(1:nZ),'reverse'); 
  % edge_effect = edge_effect - linspace(edge_effect(1), ...
  %               edge_effect(end),numel(edge_effect))';

  taper = max(0, abs(time-spike_time) - (d.dt_dx * nP));
  taper = exp(-taper.^2);

  roi = (taper > 1e-2 & taper < 0.5);
  % trend = fit(time(roi),V(roi),'poly3');
  trend = (time(roi) * [1 0] + [0 1]) \ V(roi); 
  roi = (taper > 1e-3); 
  % trend = trend(time(roi));
  trend = time(roi) * trend(1) + trend(2);
  V(roi) = V(roi) - trend;

  if any(isnan(V(:)))    
    error('NaN in computed V_extracellular')
  end
  
  wave(:,:,kk) = ( V .* taper );
% wave(:,ee) = wave(:,ee) + (V+edge_effect/k) .* taper;

  if isfield(options,'plot_wave') && ii == options.plot_wave 
    %%
    clf, hold on, lc = lines(7);
    plot(time,V-vPad,'--','color',(lc(1,:)+1)/2)    
    plot(time,V,'-','color',lc(1,:),'linewidth',1)

    plot(time,vPad,'-','Color',lc(6,:))

    if exist('C1','var')      
      plot(time,-edge_effect/dx,'-','Color',(lc(2,:)+2)/3,'linewidth',1)
      plot(time,(V+edge_effect/dx) .* taper,'k-','LineWidth',1.2)
    else plot(time,V .* taper,'k-','LineWidth',1.2)        
      plot(time(roi),trend,'-','Color',(lc(2,:)+2)/3,'linewidth',1)
    end

    plot(time,0.15 * taper,':','LineWidth',1.3,'Color',lc(5,:))
    xlim(spike_time + [-10 10])
  end

end % loop over PCA components


