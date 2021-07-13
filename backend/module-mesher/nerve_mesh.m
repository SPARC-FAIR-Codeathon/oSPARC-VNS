
function nerve_mesh(varargin) 

RUN_nerve_mesh (varargin{:})





function RUN_nerve_mesh (varargin)
% nerve_mesh(geometry, ...) builds just the mesh for a nerve and electrode
% geometry, which is saved as a EIDORS mesh in a .mat file. 
% 
% This process consists of meshing the model in GMSH and generating any 
% neccessary thin-layers as a mesh postprocessing step. 
% 
% Extracted from models.nerve_anatomy, which runs the electroanatomical 
% model implemented in EIDORS for electrical stimulation & recording.
% 
% "geometry" is a struct (or filename for a file containing a struct which
% defines the stimulating / recording electrode geometry, some meshing 
% parameters, and which fascicle set to use (see mesh.read_dat_file). 
% If "geometry" is not supplied, it's read from 
%  >> tools.cache('path','elec-geom.mat')
% 
% Some additional -options can be specified: 
% -regenerate: if set, the mesh is regenerated even if
%              tools.cache('path')/pelvic_nerve-thin.msh.mat exists. 
% 
% CDE 08-Apr-2021 : stripped out the stimulation and recording part,
% leaving just the meshing part. 

%% Parse model inputs

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

e = get_input_structure(varargin{:}); 

if any(named('-out'))
     output_name = get_('-out');     
else output_name = tools.file('get','out~/eidors-mesh (%d).mat','next'); 
end

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end

% if output_name exists but does not end in .mat: 
if ~isempty(output_name) && isempty(regexp(output_name,'\.mat$','once'))
    output_name = sprintf('%s (%s).mat', tools.file('out~/eidors-mesh'), ...
                                         output_name);
                       
    printf('\nGenerating %s', tools.file('T', output_name))
end

tools.setupEIDORS;

%% Generate requested mesh
if any(named('-fix-m')), m = old_eidors_file.model;
else m = generate_nerve_mesh(get_, named, e); 
     m = configure_array_electrodes(m); % Configure electrodes. 
end

%% Show the volume and electrode pattern

if ~isdeployed && ~any(named('-no-plot'))
    plots.view_mesh(m,'-fasc')
end
% preview_eidors_mesh(m)

printf('\nSaving %s ...', tools.file('T', output_name))
save(output_name,'-struct','m')
printf('Done!\n')

return

function m = generate_nerve_mesh(get_, named, e)

PN_mesh = tools.cache('path','nerve_mesh'); 
PN_mesh = @(ext)[PN_mesh ext]; % output file

if any(named('-regenerate')) || ~exist(PN_mesh('-thin.msh.mat'),'file')

  geo_template = dir(tools.file('~/input/array/*.geo.template','-nomake')); % Check for a subject-specific file first
  
  if any(named('-array')), geo_template = get_('-array'); 
  elseif ~isempty(geo_template)
    geo_template = [geo_template(1).folder filesep geo_template(1).name];
  elseif exist(PN_mesh('.geo.template'),'file') % try the cache 
    geo_template = PN_mesh('.geo.template');
  else error('Could not find ~/input/array/*.geo.template')
  end
  
  t = tic;
  
  make_opts = {'-usev','-output',PN_mesh('.geo'), ...
                 'virtual_thinlayer','exterior_len',e.mesh.MeshLengthMax};
  
  if any(named('-make-o'))    
    more_opts = get_('-make-o');
    if ~iscell(more_opts), more_opts = {more_opts}; end
    make_opts = [make_opts more_opts];
  end
  
  if ~exist(PN_mesh('.geo'),'file')
    tools.makeFromTemplate(geo_template, make_opts{:})
  else   warning('pnModel:existing_geo_file', ...
                 'Using existing .geo file (%s). %s', PN_mesh('.geo'), ...
                 'If this is not intended, call tools.cache(''reset'')')
  end
  
  % use gmsh.exe to convert the .geo file to a .msh file:
  if ~exist(PN_mesh('.msh'),'file')  
    path_to_gmsh = tools.configuration('gmsh'); % 'C:\Program Files\gmesh\gmsh.exe';
    system(sprintf('"%s" "%s" -3',path_to_gmsh,PN_mesh('.geo'))); 
  end

  if ~exist(PN_mesh('.msh'),'file')
    % winopen(PN_mesh('.geo'))
    error('Pipeline crashed during mesh generation')
  end
  
  % Add thinlayer and convert to EIDORS forward model (takes a minute or two)
  
  m = mesh.make_gmsh_thinLayer(PN_mesh('.geo'), e.nerve.Perineurium_mm );
  fprintf('Elapsed time is %.0f:%02.3f\n', floor(toc(t)/60), mod(toc(t),60))  
else m = load(PN_mesh('-thin.msh.mat')); 
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

%% Load the array configuration from .mat or .json file; check defaults 
function e = get_input_structure(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

is_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if nargin > 0 && exist(varargin{1},'file'), file = varargin{1}; 
elseif any(named('-geom')), file = get_('-geom'); 
elseif exist(tools.cache('path','elec-geom.mat'),'file'),
        file = tools.cache('path','elec-geom.mat'); 
elseif exist(tools.cache('path','elec-geom.json'),'file'),
        file = tools.cache('path','elec-geom.json');
else    file = tools.file('~/input/array/default.json');
    warning('pnModel:default_config_file','Using default file: %s', tools.file('T',file))
end

if isstruct(file), e = file; 
else   
  if is_ext_(file,'.mat'),      e = load(file); 
  elseif is_ext_(file,'.json'), e = tools.parse_json(file);
  elseif is_ext_(file,'.xml'),  e = tools.parseXML(file); 
  else error('unknown filetype on file "%s", expected {.mat, .json, .xml}', file)
  end

  if ~any(named('-q')), fprintf('Loading %s\n',file); end
  % CE 9-4-21: I haven't tested the .json or .xml here yet. 
  
end

%% Add some default fields 

if ~isfield(e,'mesh'), e.mesh = struct; end
if ~isfield(e.mesh,'MeshLengthMax'), e.mesh.MeshLengthMax = 0.1; end

if ~isfield(e,'nerve'), e.nerve = struct; end
if ~isfield(e.nerve,'Perineurium_mm'), e.nerve.Perineurium_mm = 1e-4; end
if isfield(e.nerve,'Perineurium_um'),  e.nerve.Perineurium_mm = 1e-3 * ...
                                       e.nerve.Perineurium_um; 
end


if ~any(named('-no-cache')), 
    save(tools.cache('path','elec-geom.mat'),'-struct','e')
end


return
