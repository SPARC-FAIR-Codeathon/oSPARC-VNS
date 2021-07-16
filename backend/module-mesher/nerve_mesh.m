
function nerve_mesh(array_file,nerve_xml,nerve_script, varargin)
% main( array.json , nerve.xml , [nerve-mesh-opts.json], [...] )
% 
% run models.nerve_mesh
%  
% main() with 0 inputs uses ~/input/array/default.json
% main( array_file ) uses the input .mat or .json file
% main( array_file, nerve_xml ) additionally uses the specified nerve
%       anatomy file (over-riding any file embedded in the config file).
% main( array_file, nerve_xml, nerve_script ) and
% main( array_file, [],        nerve_script ) uses the specified script
%       ( either .mat  or .json ) to set up the input nerve geometry.
% 
% main( ..., -flags ) exposes some additional functionality: 
% 
%   -p preview fascicle pattern to be meshed
%   -u unit test GMSH strings from input files
%   -r fast mode, do not clear cache (useful for debugging)
%   -a additional args, everything past this point is passed into
%       models.nerve_mesh as an additonal input argument
% 
% Calvin Eiber 19-April-2021

if nargin < 1 || isempty(array_file)
  fprintf('Arg 1 not set, using default array.json file\n')
  array_file = './input/demo/array.json'; 
end
if nargin < 2
  fprintf('Arg 2 not set, using default nerve.xml file\n')
  nerve_xml = './input/demo/nerve.xml'; 
end
if nargin < 3, nerve_script = ''; 
elseif strcmp(nerve_script,'demo')
  nerve_script = './input/demo/nerve-script.json'; 
end

named = @(v) strncmpi(v,varargin,length(v)); 
has_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end


if true
    disp('====================================================')
    printf('Running models.nerve_mesh %s\n', datestr(now))
    printf('{1} = %s\n{2} = %s\n{3} = %s\n',array_file, nerve_xml, nerve_script)
    disp('====================================================')
end

tools.file('root',pwd); % set 'root' to this folder

if isdeployed
  if exist(array_file,'file')
       printf('Found array description file: %s\n', array_file)
  else printf('Failed to find array description file: %s\n', array_file)
  end
end

printf('Loading %s', array_file);
if has_ext_(array_file,'.mat'), opts = load(array_file); 
elseif has_ext_(array_file,'.json'), opts = tools.parse_json(array_file);
else error('%s: unknown extension, expected .mat or .json', array_file)
end      
% can be [./array, ./mesh, ./nerve], or [array/...] internally

if ~isfield(opts,'array'), opts.array = opts; end
if isfield(opts.array,'mesh'), opts.mesh = opts.array.mesh; end

if ~isempty(nerve_script)
 
  printf('Loading %s', nerve_script);
  if has_ext_(nerve_script,'.mat'), n = load(nerve_script); 
  elseif has_ext_(nerve_script,'.json'), n = tools.parse_json(nerve_script);
  else error('%s: unknown extension, expected .mat or .json', nerve_script)
  end      

  if isfield(n,'nerve'), opts.nerve = n.nerve; else opts.nerve = n; end
  if isfield(n,'mesh'), opts.mesh = n.mesh; end
end

if ~isempty(nerve_xml), opts.nerve.file = nerve_xml;
elseif ~isfield(opts,'nerve'), opts.nerve.file = nerve_xml;
end

if any(named('-m'))
  % test mesh.insert_ ... components
  disp('GMSH code unit test')
  s = mesh.insert_gmsh_fascicles(opts);
  disp(s)
  
  a = mesh.insert_gmsh_electrodes(opts);    
  disp(a)
  return
end

%% Core code

t = tic; 
if isdeployed, disp('Progress: 3%'), end

mesh.insert_gmsh_fascicles('-setup',opts);
if ~isfield(opts.nerve,'Perineurium_um')
 try opts.nerve.Perineurium_um = mesh.get_perineurium_thickness; end
end

if any(named('-p')), plots.preview_fascicles, 
    title('click to continue or esc to cancel')
    [~,~,b] = ginput(1);
    if any(b==27), error('Cancelled'), end
end
if isdeployed, debug_path_config_for_oSPARC, end

if ~any(named('-r')), tools.cache('reset'), end

if any(named('-a')), CLI_args = varargin(find(named('-a'))+1:end); 
else CLI_args = {}; 
end

RUN_nerve_mesh('-geom',opts,CLI_args{:}, ...
               '-out',tools.file('get','./output/eidors-mesh (%d).mat','next'))

if ~any(named('-q')) || isdeployed, toc(t), end

return

%% Example unit test

function debug_path_config_for_oSPARC

fprintf('/*************************************** \n');
fprintf('  oSPARC model configuration \n');
fprintf('  pwd = "%s"\n', pwd);
fprintf('  tools.file(''~'') = "%s"\n', tools.file('~'));
fprintf('  tools.file(''out~'') = "%s"\n', tools.file('out~'));
fprintf('  tools.configuration(''noload'') = \n')
d = tools.configuration('noload'); 
disp(d)
fprintf('  tools.cache = "%s"\n', tools.cache);
fprintf(' ***************************************/ \n');






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
varargin = tools.opts_to_args(varargin,'array','mesh');
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

printf('try tools.setupEIDORS');
tools.setupEIDORS; % if ~isdeployed, tools.setupEIDORS; end

if isdeployed, disp('Progress: 6%'), end

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

if isdeployed, disp('Progress: 100%'), end


return

function m = generate_nerve_mesh(get_, named, e)

PN_mesh = tools.cache('path','nerve_mesh'); 
PN_mesh = @(ext)[PN_mesh ext]; % output file

if any(named('-regenerate')) || ~exist(PN_mesh('-thin.msh.mat'),'file')

  % Check for a subject-specific file first
  G = dir(tools.file('~/input/*.geo.template','-nomake')); 
  
  if any(named('-array')), G = get_('-array'); 
  elseif any(named('-template')), G = get_('-template');
  elseif isfield(e,'array') && isfield(e.array,'Template'), 
    G = e.array.Template;      
  elseif isfield(e,'array') && isfield(e.array,'template')
    G = e.array.template;
  elseif ~isempty(G), [~,idx] = max([G.datenum]); 
    G = [G(idx).folder filesep G(idx).name];
  elseif exist(PN_mesh('.geo.template'),'file') % try the cache 
    G = PN_mesh('.geo.template');
  else error('Could not find ~/input/*.geo.template')
  end
  
  if ~any(G == '.'), G = [G '.geo.template']; end
  if any(G == '~'), G = tools.file(G); end
  
  t = tic;
  
  make_opts = {'-usev','-output',PN_mesh('.geo'), ...
                 'virtual_thinlayer','exterior_len',e.mesh.MeshLengthMax};
  
  if any(named('-make-o'))
    more_opts = get_('-make-o');
    if ~iscell(more_opts), more_opts = {more_opts}; end
    make_opts = [make_opts more_opts];  
  end
  
  if isfield(e,'array') && isfield(e.array,'carrier')
    more_opts = tools.opts_to_args({e.array},'carrier','--s2a-keep');
    make_opts = [make_opts more_opts];  
  end
  
  if ~exist(PN_mesh('.geo'),'file')
    tools.makeFromTemplate(G, make_opts{:})
  else   warning('pnModel:existing_geo_file', ...
                 'Using existing .geo file (%s). %s', PN_mesh('.geo'), ...
                 'If this is not intended, call tools.cache(''reset'')')
  end

  if isdeployed, disp('Progress: 10%'), end
  path_to_gmsh = tools.configuration('gmsh'); % 'C:\Program Files\gmesh\gmsh.exe';

  if any(named('-preview-geo')), 
      m = PN_mesh('.geo'); 
     system(sprintf('"%s" "%s" &',path_to_gmsh,PN_mesh('.geo')));
     return
  end
  
  % use gmsh.exe to convert the .geo file to a .msh file:
  if ~exist(PN_mesh('.msh'),'file')  
    system(sprintf('"%s" "%s" -3',path_to_gmsh,PN_mesh('.geo'))); 
  end

  if isdeployed, disp('Progress: 30%'), end
  
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
