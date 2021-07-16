
function [gmsh,config] = insert_gmsh_fascicles(varargin)
% gmsh = insert_gmsh_fascicles inserts fascicles (from a MBF .xml file) 
%   into a GMESH .geo file during the compilation process
%   (using make_from_template, usually), such that the fascicles align
%   to the X-axis in GMESH (this usually comes out as the Z-axis later)
%   example: 
% @[insert_gmsh_fascicles] // Execute spline script using SETUP
% 
% Syntax: -setup           : adds persistent options to static gmsh call
%                           (this is used to customise the behaviour of 
%                            a makeFromTemplate call) 
%         -anat [filename] : select source file, default is: 
%                            input~/nerve/sub-131_sam-0_fascicleTracing.xml
%         -anat [struct]   : use preloaded anatomy source structure
%         -conf [filename] : select configuration file, default is: 
%                             tools.cache('nerve-config.mat')
%         -conf [struct]   : use preloaded configuration structure
% 
% The additional following sets of options are defined:
%  Fascicles can be re-arranged to make test patterns using the following: 
%   -pIndex [1..n]        : select fascicle indices from source file
%   -pLocation [1..n x 2] : set fascicle XY location
%   -pRotate [1..n]       : set rotation (in degrees) for each fascicle
%                           around fascicle centroid
% 
%  The overall resulting pattern of fascicles can be moved, shifted and
%    scaled using the following: 
%   -xMove [x y]; -xMove [y]
%   -xScale [xy]; -xScale [x y]; -xScale [2x2 matrix]
%   -xRotate [deg];
%  These translations, scalings, and rotations are applied around the
%    centroid of the overall fascicle pattern and are applied in the order 
%    [scale, rotate, move]
% 
% By default, the fascicles are extruded from -z to +z with z defined using
%   -z [range]. More complex fascicle trajectories can be specified using
%   -path [xyz], where the xyz path of the fascicle is specified in mm. 
% 
% Input units are assumed um and the output units are assumed mm. 
%  for input units of mm, use -units-mm
%  for output units of um, use -units-um
% 
% v0.5, Calvin Eiber 9 April 2021 major refactor for oSPARC 
% v0.4, Calvin Eiber renamed from insert_gmsh_splines
% v0.3, Calvin Eiber 1 May 2020 - major re-factor, implemented this syntax.
%                     warning: not all combinations of options implemented. 
% 
persistent persistent_argin info_result

if any(strcmpi(varargin,'-setup'))
  persistent_argin = varargin(~strcmpi(varargin,'-setup'));
  info_result = []; 
elseif iscell(persistent_argin)
  varargin = [persistent_argin varargin]; 
end 

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1}; %#ok<NASGU>

if any(named('-setup')) && nargout == 0, return, end
if any(named('-clear')), info_result = []; return, end

if ~isempty(info_result) && (any(named('-info')) || nargout > 1)
  gmsh = info_result; return % Return the processed files and configuration    
end

[anat,config] = parse_input_arguments(varargin{:}); 

if any(named('-setup')) % get confirmation info and exit 
  if nargout > 1, gmsh = anat; 
  elseif nargout > 0, gmsh = src_file; 
  end, return
end

if any(cellfun(@(x)isfield(config,x),{'pIndex','pLocation'}))
	anat = construct_fascicle_pattern(anat,config,varargin{:});
end

if any(cellfun(@(x)isfield(config,x),{'xMove','xScale','xRotate'}))    
	anat = apply_geometric_transforms(anat,config,varargin{:});
end

if any(named('-info')) || nargout > 1 
  info_result = anat; 
  
  if isfield(info_result,'splines')
    if any(named('-all')) % cell output structure
      for f = fieldnames(anat.splines)'
        info_result.(f{1}) = {anat.splines.(f{1})};
      end
    else % scalar output structure (fascicles only)
      for f = fieldnames(anat.splines)'
        info_result.(f{1}) = anat.splines(1).(f{1});
      end        
    end
  end
  
  if any(named('-check-units')), anat = convert_units(anat,config,named); end
  gmsh = anat; return % Return the processed files and configuration    
end

anat = convert_units(anat,config,named);

%% INITIALISE code to add to gmsh document
gmsh = sprintf('\n// Fascicles (from %s)\n', anat.filename);

this = struct; 
this.resol = config.MeshLengthFascicle; 

w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:});
gmsh = w_(gmsh,'lc = %0.6f; // mesh resolution',   this.resol   );
gmsh = w_(gmsh,'f_id = newv;\n' );

insert_(anat.splines,'reset'); 

nF = size(anat.splines(1).coeffs,3);

%% FOR EACH compartment 

for ii = numel(anat.splines):-1:1
    
  this.DoSurface = strcmpi(anat.splines(ii).type,'Fascicle');
  gmsh = insert_(gmsh,named,anat.splines(ii),this);
  
end

return

%% 
function [source_anat,config] = parse_input_arguments(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end

%% Parse source file (-anatomy)

if any(named('-anat')), source_file = get_('-anat'); 
elseif any(named('-file')), source_file = get_('-file'); 
elseif nargin > 0 && isstruct(varargin{1}), source_file = varargin{1};     
elseif nargin > 0 && exist(varargin{1},'file'), source_file = varargin{1};
elseif nargin > 0 && any(varargin{1} == '~'), source_file = tools.file(varargin{1});
else      
  error('no source file for fascicle pattern')
end

has_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if isstruct(source_file) % precomputed 
  source_anat = source_file; 
  
  % input may be the configuration structure, use that if relevent. 
  if isfield(source_anat,'nerve')
      config_file = source_anat; 
      if isfield(source_anat.nerve,'file')
           source_file = {'-anat',source_anat.nerve.file};
        if any(source_file{2} == '~'), source_file{2} = tools.file(source_file{2}); end
      else source_file = {'-anat',config_file.nerve};
      end
      [source_anat,config] = parse_input_arguments(source_file{:}, ...
                                                   '-conf',config_file, ...
                                                      varargin{:});
      return
  end % if (config structure)
  
elseif has_ext_(source_file,'.mat') % from mat file
  printf('Loading %s', source_file);
  source_anat = load(source_file);  
  source_anat.filename = source_file;
elseif has_ext_(source_file,'.xml') % from XML
  printf('Loading %s', source_file);
  source_anat = tools.parse_xml(source_file);  
  source_anat.filename = source_file;
  source_anat.Type = 'MBF-XML';  
  
elseif has_ext_(source_file,'.json') % from XML
    
  printf('Loading %s', source_file);
  source_anat = tools.parse_json(source_file);
  source_anat.filename = source_file;  
  source_anat.Type = 'JSON';
  
  % input may be the configuration structure, use that if relevent. 
  if isfield(source_anat,'nerve') && isfield(source_anat.nerve,'file')
      config_file = source_anat; 
      source_file = source_anat.nerve.file;         
      [source_anat,config] = parse_input_arguments('-anat',source_file,'-conf',config_file,varargin{:});
      return
  end % if (config structure)
  
elseif has_ext_(source_file,'.splines.dat')  
  printf('Loading %s', source_file);
  
  index = []; 
  
  if any(named('-index')), index = get_('-index');
  elseif any(named('-con')), config_file = get_('-con'); 
    if exist(config_file,'file') && has_ext_(config_file,'.mat')
      config_file = load(config_file);
    end
      
    if isstruct(config_file)
      if isfield(config_file,'index'), index = config_file.index;
      elseif isfield(config_file,'SplineIndex'), index = config_file.SplineIndex;      
      end
    else clear config_file % because exist(var) used below
    end
  else index = []; 
  end

  if isempty(index)
     warning('insert_gmsh_splines:defaultIndex', '-index not set, using default (#1)')
            index = 1;
  end
  
  d = mesh.read_dat_file(source_file,'index',index);
 
  source_anat.index = index; 
  source_anat.filename = source_file; 
  source_anat.splines.coeffs = d.coeffs; 
  source_anat.splines.outline = d.outline; 
  source_anat.splines.type = {'Fascicle'}; 
  source_anat.splines.info = sprintf('Imported from %s', source_file);
  
  if isfield(d,'info') % && ~iscell(d.info)
    source_anat.info = d.info; % strsplit(d.info,newline); 
  end
else
  error('unknown filetype on file "%s", expected {.mat, .xml}', source_file)
end

%% Parse config file (-config)
if exist('config_file','var') % passed in as argument or already loaded
   if isdeployed, printf('Using input configuration structure'); end
elseif any(named('-con')), config_file = get_('-con'); 
elseif nargin > 2 && ischar(varargin{1}) && ischar(varargin{2}) && ...
                exist(varargin{1},'file') && exist(varargin{2},'file')
     config_file = varargin{2};
elseif exist(tools.cache('path','nerve-config.mat'),'file')
     config_file = tools.cache('path','nerve-config.mat'); 
else config_file = struct; 
  if nargin < 2, printf('[Warning] no configuration file found.'); end
end

if isstruct(config_file) 
    
  config = config_file; 
  if isfield(config,'nerve')
    if isfield(config,'mesh') % union of 'nerve' and 'mesh' properties
        c = config.mesh; 
        for f = fieldnames(config.nerve)'
            c.(f{1}) = config.nerve.(f{1});
        end
         config = c; 
    else config = config.nerve; 
    end
  end
elseif isempty(config_file), config = struct; % empty
elseif has_ext_(config_file,'.mat') % from mat file
  printf('Loading %s', config_file);
  config = load(config_file);  
  config.filename = config_file;
elseif has_ext_(config_file,'.xml') % from XML
  printf('Loading %s', config_file);
  config = tools.parse_xml(config_file);
  config.filename = config_file;
  warning untested
elseif has_ext_(config_file,'.json') % from JSON
  printf('Loading %s', config_file);
  config = tools.parse_json(config_file);
  if isfield(config,'nerve'), config = config.nerve; end
  config.filename = config_file;
  warning untested
else
  error('unknown filetype on file "%s", expected {.mat, .json}', file)
end

% if any(named('-i')),  config.index = get_('-i');  end
if any(named('-xM')), config.xMove = get_('-xM');   end 
if any(named('-xS')), config.xScale = get_('-xS');  end
if any(named('-xR')), config.xRotate = get_('-xR'); end

if any(named('-z')),  config.zRange = get_('-z');   
elseif any(named('-domain')), config.zRange = get_('-domain');
end

if ~isfield(config,'zRange'), config.zRange = 6; 
elseif numel(config.zRange) > 1, config.zRange = config.zRange(1); 
end

if any(named('-path')), config.FascicleTrajectory = get_('-path'); 
elseif any(named('-curve')), config.FascicleTrajectory = get_('-curve'); 
end

if any(named('-pI')), config.pIndex = get_('-pI');    end
if any(named('-pL')), config.pLocation = get_('-pL'); end
if any(named('-pR')), config.pRotation = get_('-pR');   end


if any(named('-meshList')),   config.meshList = get_('-meshList'); 
elseif any(named('-mesh-l')), config.meshList = get_('-mesh-l'); 
elseif any(named('-do')),     config.meshList = get_('-do');
elseif ~isfield(config,'meshList'), config.meshList = {'Fascicle'}; 
end

if any(named('-lc')), config.MeshLengthFascicle = get_('-lc');
elseif any(named('-res')), config.MeshLengthFascicle = get_('-res');
elseif any(named('-MeshLen')), config.MeshLengthFascicle = get_('-MeshLen');
elseif ~isfield(config,'MeshLengthFascicle')
    config.MeshLengthFascicle = 0.02; 
end

if ~isfield(source_anat,'splines') % already done
    source_anat = generate_splines(get_,named,source_anat, config); 
end

if isfield(config,'FascicleTrajectory')
    
    if isfield(source_anat,'trajectory')
      warning('insert_gmsh_splines:overrideFascicleTrajectory', ...
              '%s had trajectories but I''m using the ones %s', ...
              source_anat.filename, 'in the configuration struct.')
    end
  if iscell(config.FascicleTrajectory)
    source_anat.trajectory = config.FascicleTrajectory;       
  else
    nF = size(source_anat.splines.coeffs,3); 
    source_anat.trajectory = cell(1,nF);
    source_anat.trajectory(:) = {config.FascicleTrajectory};
    if isfield(config,'FascicleTrajectoryIntercept')
      for ff = 1:nF
        source_anat.trajectory{ff} = source_anat.trajectory{ff} - ...
                                     source_anat.trajectory{ff}(1,:) + ...
                                     config.FascicleTrajectoryIntercept(ff,:); 
      end
    end
  end
end

return



%% Convert MBF XML to splines 
function xml = generate_splines(get_,named,xml,config)

if isfield(xml,'splines'), return, end % already done
if ~isfield(xml,'Children'), % not XML
  if isfield(xml,'coeffs'),  xml.splines.coeffs  = xml.coeffs;  end
  if isfield(xml,'outline'), xml.splines.outline = xml.outline; end
  if isfield(xml,'type'),    xml.splines.type = xml.type;
  else                       xml.splines.type = {'Fascicle'}; 
  end
  return
end % not XML

%% x.Children().name representation
the = @(x,n) x.Children(strncmpi({x.Children.Name},n,length(n)));
attr_ = @(x,n) x.Attributes(strcmpi(n,{x.Attributes.Name})).Value;
trim_ = @(x) x.Children(~contains({x.Children.Name},'#text'));

s = struct; % splines accumulator

if any(named('-nc')), nCoe = get_('-nc'); else nCoe = 25; end
if any(named('-np')), nPxy = get_('-np'); else nPxy = 127; end
if any(named('-nt')), nSwe = get_('-nt'); else nSwe = 255; end

if any(named('-invert-y')), ys = 1; else ys = -1; end % fix MDF XML -Y

xml.Children = trim_(xml);
contours = the(xml,'contour');

% DEV_demo = mesh.read_dat_file; 
% DEV = evalin('base','DEV_demo'); 
%%

s.coeffs = {};
s.outline = {};
s.type = {}; 
s.info = sprintf('Generated from %s', xml.filename); 

for ii = 1:numel(contours)
    
    this = contours(ii);     
    type = regexprep(attr_(this,'name'),'-\d+','');
    if ischar(type), type = strtrim(type); end
    if ~ismember(type,s.type)
         s.type = [s.type {type}];
         tid = numel(s.type);
    else tid = find(strcmpi(s.type,type)); 
    end

    points = the(this,'point'); 
    
    x = arrayfun(@(p) str2double(attr_(p,'x')), points);
    y = arrayfun(@(p) str2double(attr_(p,'y')), points) * ys; 
    l = cumsum(sqrt((x-x([2:end 1])).^2 + (y-y([2:end 1])).^2));
    ok = diff([0 l]) > eps; x = x(ok); y = y(ok); l = l(ok); % remove duplicate points
    
    if numel(s.coeffs) < tid
        s.coeffs{tid} = []; 
        s.outline{tid} = [];
    end
    
    f = size(s.coeffs{tid},3) + ~isempty(s.coeffs{tid}); 
    
    x(end+1) = x(1); y(end+1) = y(1);  %#ok<AGROW>
    l(end+1) = 2*l(end)-l(end-1);      %#ok<AGROW>    
    
    s.coeffs{tid}(1,:,f) = interp1(l,x,linspace(l(1),l(end),nCoe));
    s.coeffs{tid}(2,:,f) = interp1(l,y,linspace(l(1),l(end),nCoe));
    s.outline{tid}(:,1,f) = interp1(l,x,linspace(l(1),l(end),nPxy));
    s.outline{tid}(:,2,f) = interp1(l,y,linspace(l(1),l(end),nPxy));
    
    extrude_trajectory = [the(this,'path') the(this,'sweep')]; % synonyms
    
    if ~isempty(extrude_trajectory)

      x = arrayfun(@(p) str2double(attr_(p,'x')), extrude_trajectory);
      y = arrayfun(@(p) str2double(attr_(p,'y')), extrude_trajectory) * ys; 
      l = cumsum(sqrt((x-x([2:end 1])).^2 + (y-y([2:end 1])).^2));
      ok = diff([0 l]) > eps; x = x(ok); y = y(ok); l = l(ok); % remove duplicate points

      this_path(:,1) = interp1(l,x,linspace(l(1),l(end),nSwe));
      this_path(:,2) = interp1(l,y,linspace(l(1),l(end),nSwe));
      if ~isfield(s,'trajectory'), s.trajectory = {}; end      
      
      if numel(s.trajectory) >= f && ~isempty(s.trajectory{f})
        warning('insert_gmsh_splines:multipleFascicleTrajectory', ...
                '%s %d had trajectories (<path x=... \\> entries), %s %s.', ...
                'Multiple objects associated with fascicle', f,  ...
                'ignoring the one associated with', type)
      end
      
      s.trajectory{f} = this_path;
    end
end

%% Put fascicles first in sequence

if numel(s.type) == 1
  s.type = s.type{1}; 
  s.coeffs = s.coeffs{1};
  s.outline = s.outline{1};
else
  for keyword = {'FascicleInterior','Fascicle'}

    k = strcmpi(s.type,keyword{1});
    if ~any(k), continue, end 
    k = [find(k) find(~k)];
    
    s.coeffs = s.coeffs(k);
    s.outline = s.outline(k);
    s.type = s.type(k); 
    break
  end
end


if numel(config.meshList) == 1 && ...
   strcmp(config.meshList{1},'Fascicle') % Default scenario
  s.type = 'Fascicle';
  s.coeffs = s.coeffs{1};
  s.outline = s.outline{1};
else
  type_index = zeros(size(config.meshList));
  for obj = 1:numel(type_index)
    sel = strcmp(s.type,config.meshList{obj});
    if any(sel), type_index(obj) = find(sel,1); continue, end
    types = lower(s.type); 
    count = cellfun(@(x) size(x,3), s.coeffs);
    
    switch config.meshList{obj}
      case {'Fascicle'}          
        sel = contains(types,'fascicle');
        if ~any(sel), sel = contains(types,'perineur'); end
        if ~any(sel), sel = contains(types,'endoneur'); end
        if any(sel & ~contains(types,'outer'))
          sel = sel & ~contains(types,'outer'); 
        end
        if sum(sel) > 1, sel = sel & count == max(count(sel)); end
        type_index(obj) = find(sel,1);
        
      case {'Epineurium'}

        sel = contains(types,'epineur');
        if ~any(sel), sel = ~contains(types,' of '); end
        if sum(sel) > 1, sel = sel & count == min(count(sel)); end
        type_index(obj) = find(sel,1);
          
      otherwise 
        error('Unknown %s "%s". Did you mean "%s" (%s)?', ...
              'meshList entry', config.meshList{obj},'Fascicle', ...
              'these are case-sensitive');
    end
        
     if type_index(obj) == 0
         tlist = sprintf(', ''%s''', s.type{:});
         error('"%s" did not match any contours from %s: {%s}', ...
                config.meshList{obj}, xml.filename, tlist(3:end))
     end
      
  end
    
  s = struct('type',config.meshList, ...
             'coeffs',s.coeffs(type_index), ...
             'outline',s.outline(type_index));
end

xml.splines = s; 

return


%% Build artificial fascicle-patterns using the -xy command
function anat = construct_fascicle_pattern(anat,config,varargin)

named = @(v) strncmpi(v,varargin,length(v)); 

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end

if isfield(anat,'source'),
     s = anat.source; 
else s = anat.splines; 
     anat.source = s; 
end, s0 = s; 

if numel(s) > 1, error('TODO multiple tissues'), end

%%
f_index = config.pIndex; 
if isfield(config,'pLocation'), f_xy = config.pLocation; 
else f_xy = reshape(0*f_index,[],1);
end
f_angle = 0; 

if iscell(f_xy), f_xy = cat(1,f_xy{:}); end
if iscell(f_xy), f_xy = reshape(cat(1,f_xy{:}),size(f_xy,1),[]); end
if size(f_xy,2) == 1, f_xy = f_xy(:) * [0 1]; end

if isempty(f_index), return, end
if isfield(config,'pRotation'), f_angle = config.pRotation; end

nF = max(numel(f_index),size(f_xy,1)); 

if numel(f_index) < size(f_xy,1), f_index = repmat(f_index(1),[1 nF]); end
if numel(f_angle) < nF, f_angle = repmat(f_angle(1),[1 nF]); end

if iscell(f_index), f_index = [f_index{:}]; end

printf('Building synthetic fascicle pattern (%d fascicles)', nF);

%%

% Start building input structure 

if ~isfield(s,'info') && isfield(s,'filename'), s.info = s.filename;
elseif ~isfield(s,'info'), s.info = 'MATLAB defined geometry'; 
end

s.info = sprintf('XY pattern from %s', s.info);
s.source = s.coeffs; 
s.fascicle_id = f_index;
s.fascicle_xy = f_xy; 

s.coeffs = []; 
s.outline = []; 

if iscell(s0.outline), 
  s0.outline = s0.outline{1}; 
  s0.coeffs  = s0.coeffs{1};     
  if any(named('-do')), error TODO_for_perineurium_etc, end
end


for ii = 1:numel(f_index), ff = f_index(ii);
    
    printf('Fascicle %d: input %d, xy = [%g %g], a = %g deg', ii, ff, f_xy(ii,:), f_angle(ii))
    
    xy0 = median(s0.outline(:,:,ff)) ; % outline already in mm 
    
    a = deg2rad(f_angle(ii)); 
    
    Q = [cos(a) sin(a); -sin(a) cos(a)];
    
    s.coeffs(:,:,ii) = Q'*(s0.coeffs(:,:,ff) - xy0') + f_xy(ii,:)';    
    s.outline(:,:,ii) = (s0.outline(:,:,ff) - xy0)*Q + f_xy(ii,:);
    
    if any(s.outline(:,2,ii) < 0) && ~isfield(anat,'trajectory')
        warning('insert_gmsh_splines:NegativeXY', sprintf(...
                'Negative Y in fascicle %d (#%d, <%0.3f %0.3f>)', ...
                ii, ff, f_xy(ii,:))) %#ok<SPWRN>
    end
end

if isfield(anat,'trajectory'), s.trajectory = anat.trajectory(f_index); 
end

anat.splines = s; 

%% apply move, scale, rotate to fascicle pattern
function anat = apply_geometric_transforms(anat,config,varargin)

s = anat.splines;
named = @(v) strncmpi(v,varargin,length(v)); 

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end

if numel(s) > 1
  new_s = [];
  for ii = 1:numel(s)      
    anat.splines = s(ii);
    this = apply_geometric_transforms(anat,config, varargin{:});    
    if isempty(new_s), new_s = this.splines;
    else new_s(ii) = this.splines; %#ok<AGROW>
    end
  end
  anat.splines = new_s;
  return
end

% From here, we can safely assume that s is a scalar structure
if iscell(s.type), s.type = s.type{1}; end

nF = size(s.coeffs,3);
xy0 = mean(s.coeffs(:,:),2)'; 

if any(named('-debug-X'))
  cla, hold on
  for ff = 1:nF      
    plot(s.outline(:,1,ff),s.outline(:,2,ff),':','LineWidth',1.2);
  end
end

if isfield(config,'xRotate')
  
  while iscell(config.xRotate), config.xRotate = config.xRotate{1}; end
  a = deg2rad(config.xRotate(1));
  Q = [cos(a) sin(a); -sin(a) cos(a)];

  printf('Rotating %s pattern by %g deg', lower(s.type), config.xRotate(1))
  
  for ff = 1:nF
    s.coeffs(:,:,ff) = Q'*(s.coeffs(:,:,ff) - xy0') + xy0' ;
    s.outline(:,:,ff) = (s.outline(:,:,ff) - xy0)*Q + xy0  ;
  end
end

if isfield(config,'xScale')
    
  while iscell(config.xScale), config.xScale = [config.xScale{:}]; end
  switch numel(config.xScale)
      
    case 1, Q = eye(2) * config.xScale; 
    case 2, Q = diag(config.xScale);
    case 4, Q = reshape(config.xScale,2,2);
    otherwise error('xScale: expected 1, 2, or 4-element vector or matrix')
  end
  
  printf('Scaling %s pattern by [%g %g; %g %g]',lower(s.type), Q')
  
  for ff = 1:nF
    s.coeffs(:,:,ff) = Q'*(s.coeffs(:,:,ff) - xy0') + xy0' ;
    s.outline(:,:,ff) = (s.outline(:,:,ff) - xy0)*Q + xy0  ;
  end
  
end

if isfield(config,'xMove')
    
  while iscell(config.xMove), config.xMove = [config.xMove{:}]; end

  dxy = reshape(config.xMove,1,[]); 
  if numel(dxy) == 1, dxy = [0 dxy]; end
    
  printf('Displacing %s pattern by [%g %g]',lower(s.type), dxy)
  
  for ff = 1:nF
    s.coeffs(:,:,ff) = (s.coeffs(:,:,ff) + dxy') ;
    s.outline(:,:,ff) = (s.outline(:,:,ff) + dxy);
  end
end

if any(named('-debug-X'))
  C = lines(nF); 
  for ff = 1:nF
    plot(s.outline(:,1,ff),s.outline(:,2,ff),'-','LineWidth',1.2,'Color',C(ff,:));
  end
  axis image, tools.tidyPlot
end

anat.splines = s; 

return

%% convert units from um to mm
function anat = convert_units(anat,config,named)

if ~isfield(anat,'splines') && isfield(anat,'coeffs')
  if ~isfield(anat,'type'), anat.type = 'Fascicle'; end
  anat.splines = anat;
end

if ~isfield(anat,'filename')
  if isfield(anat,'info'), anat.filename = anat.info(1:min(30,end));
  else anat.filename= ''; 
  end
end

if any(named('-units-um')), factor = 1; % output units um
else factor = 1000; % output mm, input um
end

if any(named('-units-mm')), factor = factor / 1000; end % input mm

for ss = 1:numel(anat.splines) % For each object class in spline file 
  anat.splines(ss).coeffs = anat.splines(ss).coeffs / factor; 
  anat.splines(ss).outline = anat.splines(ss).outline / factor;
end
[anat.splines.z] = deal(config.zRange); 





%% Write the gmsh code for each compartment 
function gmsh = insert_(gmsh, named, this, object)

persistent cross_section % Keep track of object indices, offset from vol_id
if ischar(named) && strcmp(named,'reset'), cross_section = {}; return, end
if isempty(cross_section), cross_section = {}; end

%% Locate source data

xy = this.coeffs;
w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:});

%%

for ii = 1:size(this.coeffs,3)
    
    if isempty(xy), 
      warning('insert_gmsh_splines:Empty','Possible empty %s set', ...
                                        lower(this.type))
      break
    end
  
    is_member = zeros(size(xy,2),1);
    for pp = numel(cross_section):-1:1        
        sel = in_loop(xy(:, is_member == 0, ii), cross_section{pp});
        is_member(sel) = pp; 
    end
    if false && any(is_member)
        %% check in_loop result
        clf, hold on
        for pp = 1:numel(cross_section)
            if any(is_member == pp), ls = '-'; else ls = '--'; end
            plot(cross_section{pp}(1,:),cross_section{pp}(2,:),ls)
        end
        plot(xy(1,:,ii),xy(2,:,ii),'-k','LineWidth',1.3)
        axis equal, tools.tidyPlot, grid on
    end
    cross_section = [cross_section {xy(:,:,ii)}]; %#ok<AGROW>
    is_member(is_member>0) = is_member(is_member>0); 
    
    % Check if fascicles are specified clockwise or CCW, needed to
    % calculate thin layers correctly    
    if ~any(named('-skip-ccw-check'))    
      is_ccw = ( xy(1,3:end,ii)  -  xy(1,2:end-1,ii) ) .* ...
               ( xy(2,2:end-1,ii) - xy(2,1:end-2,ii) ) -  ...
               ( xy(1,2:end-1,ii) - xy(1,1:end-2,ii) ) .* ...
               ( xy(2,3:end,ii)  -  xy(2,2:end-1,ii) );
    
      if mean(sign(is_ccw)) > 0, 
        warning('insert_gmsh_splines:WindingOrder','Flipping the order of Fascicle %d', ii)
        xy(:,:,ii) = xy(:,end:-1:1,ii);
      end
    end    
    
    vol_tag = 'su[1]';
    
    if isfield(this,'trajectory')
      
      if iscell(this.trajectory)
           this.e_this = this.trajectory{ii};
      else this.e_this = this.trajectory + [0 mean(xy(:,2:end,ii)')]; %#ok<UDIM>
      end      
      
      sstr = sweep_along_curve(this,this.type,xy(:,:,ii),ii); 
      gmsh = w_(gmsh,'%s\n',sstr);
      
      vol_tag = 'su[0]';
    elseif isfield(object,'compound_f')      && ...
            ischar(object.compound_f)        && ...
            strcmp(object.compound_f,'auto') 
          
      if test_fascicle_closeness(this, object, ii)
        sstr = build_compound_structure(this,object,this.type,xy,ii); 
        gmsh = w_(gmsh,'%s\n',sstr);
      else
        sstr = default_linear_extrude(this, this.type, xy, ii); 
        gmsh = w_(gmsh,'%s\n',sstr);
      end

    elseif isfield(object,'compound_f')      
      if any([object.compound_f{:}] == ii)
           sstr = build_compound_structure(this,object,this.type,xy,ii); 
      else sstr = default_linear_extrude(this, this.type, xy, ii); 
      end
      gmsh = w_(gmsh,'%s\n',sstr);
      
    else % Default : linear extrusion
      sstr = default_linear_extrude(this, this.type, xy, ii); 
      gmsh = w_(gmsh,'%s\n',sstr);
    end
    
    myID = numel(cross_section);
    gmsh = w_(gmsh,'ne%02d = %s;', myID, vol_tag);
    gmsh = w_(gmsh,['BooleanDifference{ Volume{vol_id}; Delete; }' ... 
                                     '{ Volume{%s}; }'],vol_tag); 
    is_member = unique(is_member(is_member>0));
    % if isfield(object,'intersections') select and append to ismember 
    for pp = 1:numel(is_member)            
        gmsh = w_(gmsh,['BooleanDifference{ Volume{ne%02d}; Delete; }' ... 
                                         '{ Volume{%s}; }'], ... 
                                            is_member(pp), vol_tag); 
    end
    if object.DoSurface
        gmsh = w_(gmsh,'Physical Surface("P_%s%d") = {su[2]};', this.type, ii);
    end
    
    gmsh = w_(gmsh,'Physical Volume("%s%d") = {%s};\n', this.type, ii, vol_tag );
end

return

%% Extrusion code

function gmsh = sweep_along_curve(splines,name,XY,id)

%% Extrude { extrude-list } Using Wire { expression-list }
norm_ = @(v) v./sqrt(sum(v.^2,2)); % normal vector

xyP = splines.e_this;              % PATH xy
xyF = [ 0*XY(1,:)'  fliplr(XY') ]; % OUTLINE xy
xyF = xyF - mean(xyF(2:end,:));    

dir = norm_(xyP(2,:)-xyP(1,:)); 
dir = [dir; norm_(cross(dir,[0 0 -1])); norm_(cross(dir,[0 1 0]))];

xyR = xyF * dir + xyP(1,:); 

if 0 
  %% Debug visualisation
  clf, hold on %#ok<UNRCH>
  plot3(xyF(:,1), xyF(:,2),xyF(:,3),'clipping','off'), grid on
  plot3(xyP(:,1), xyP(:,2),xyP(:,3),'clipping','off')
  plot3(xyR(:,1), xyR(:,2),xyR(:,3),'clipping','off'), axis equal

  a = xyP(1,:); w = 0.8; 
  axis(a([1 1 2 2 3 3]) + [-w w -w w -w w])
end

%% Add geometry to GMSH string

gmsh = ''; 
w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:});

sstr = ''; 
gmsh = w_(gmsh,'// %s_%02d\np = newp;',name,id);
for pp = 1:size(xyR,1)-1 % last point is duplicate
    gmsh = w_(gmsh, 'Point(p+%d) = {%0.24f,%0.24f,%0.24f, lc} ;', ...
                              pp-1, xyR(pp,:));
    sstr = sprintf('p+%d,%s',pp-1,sstr);
end
sstr = strrep(sstr,'p+0,','p');

gmsh = w_(gmsh,'Printf("Spline %s_%02d : %%g - %%g", p, p+%d) ;',name,id,pp-1);
gmsh = w_(gmsh,'\ne = newreg;\nSpline(e) = {p,%s} ;',sstr);
gmsh = w_(gmsh,'Curve Loop(e+1) = {e} ;');
gmsh = w_(gmsh,'Plane Surface(e+2) = e+1 ;');

% And again for extrusion path
sstr = '';
gmsh = w_(gmsh,'// %s_%02d_PATH\np = newp;',name,id);
for pp = 1:size(xyP,1)
    gmsh = w_(gmsh, 'Point(p+%d) = {%0.6f,%0.6f,%0.6f, lc} ;',pp-1,xyP(pp,:));
    sstr = sprintf('%s,p+%d',sstr,pp-1);
end
sstr = strrep(sstr,',p+0,','p,'); % also remove leading ","

gmsh = w_(gmsh,'Printf("Wire Path %s_%02d : %%g - %%g", p, p+%d) ;',name,id,pp-1);
gmsh = w_(gmsh,'\nep = newreg;\nSpline(ep) = {%s} ;',sstr);
gmsh = w_(gmsh,'Wire(ep+1) = {ep} ;');
            
gmsh = w_(gmsh,'su[] = Extrude{ Surface{e+2}; } Using Wire {ep+1} ;');
gmsh = w_(gmsh,'Printf("Extrude SU = %%g",su[0]) ;');

return

function gmsh = default_linear_extrude(splines, name, XY, xy_index) 

  gmsh = '';
  w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:});

  sstr = '';     
  gmsh = w_(gmsh,'// %s_%02d\np = newp;',name,xy_index);
  for pp = 1:size(XY,2)-1 % last point is duplicate 
      gmsh = w_(gmsh, 'Point(p+%d) = {%0.1f,%0.6f,%0.6f, lc} ;', ...
                                pp-1, -splines.z, XY(2,pp,xy_index), ...
                                                  XY(1,pp,xy_index));
      sstr = sprintf('p+%d,%s',pp-1,sstr);
  end
  sstr = strrep(sstr,'p+0,','p');

  gmsh = w_(gmsh,'Printf("Spline %s_%02d : %%g - %%g", p, p+%d) ;',name,xy_index,pp-1);
  gmsh = w_(gmsh,'\ne = newreg;\nSpline(e) = {p,%s} ;',sstr);
  gmsh = w_(gmsh,'Curve Loop(e+1) = {e} ;');
  gmsh = w_(gmsh,'Plane Surface(e+2) = e+1 ;');
  gmsh = w_(gmsh,'su[] = Extrude{%0.1f,0,0}{ Surface{e+2}; };', 2*splines.z);

function do_merge = test_fascicle_closeness(splines, info, f_id)

nF = size(splines.coeffs,3);
min_dist = nan(nF,1);

for f2 = 1:nF
  if f2 == f_id, continue, end
  d = (splines.outline(:,1,f_id)-splines.outline(:,1,f2)').^2 + ...
      (splines.outline(:,2,f_id)-splines.outline(:,2,f2)').^2 ;
    min_dist(f2) = sqrt(min(d(:))); 
end
 
do_merge = any(min_dist < info.resol);

return

function gmsh = build_compound_structure(splines,info,name,~,xy_id) 

persistent compound_structure
if ischar(splines) && strcmp(splines,'reset')
  compound_structure = []; return, 
end

if ~isempty(compound_structure)
  sel = cellfun(@(ids) any(ids == xy_id), {compound_structure.ids}); 
  if any(sel), this = compound_structure(sel); 
  else this = [];
  end
else this = []; 
end

debug_plot = true; 


if isempty(this)
  
  splines.coeffs(:,end,:) = []; % remove duplicate endpoint

  this = struct;
  this.obj = name;

  nF = size(splines.coeffs,3);  

  if isfield(info,'compound_f')
    sel = cellfun(@(x) any(x == xy_id), info.compound_f);
    this.ids = reshape([info.compound_f{sel}],1,[]);
  else % detemine automatically
    
    min_dist = nan(nF);
    for f1 = 1:nF
      for f2 = 1:nF
        if f2 == f1, continue, end
        d = (splines.outline(:,1,f1)-splines.outline(:,1,f2)').^2 + ...
            (splines.outline(:,2,f1)-splines.outline(:,2,f2)').^2 ;
        min_dist(f1,f2) = sqrt(min(d(:))); 
      end
    end

    check = false(nF,2); check(xy_id,:) = true;

    while any(check(:,1))

      check(:,1) = check(:,1) | any(min_dist(:,check(:,1)) < 0.2*info.resol,2);
      check(:,1) = check(:,1) & ~check(:,2); 
      check(:,2) = check(:,1) | check(:,2);
    end
  
    this.ids = reshape(find(check(:,2)),1,[]);
  end
  
  %%
  this.splines = {};  
  this.p_corner = []; 
  this.s_corner = zeros(0,2); 
  this.s_order = cell(nF,1); 
  this.s_gmsh_ids = {};
  this.p_gmsh_ids = {};  
  this.s_written = false;
  this.c_loops = {};
  
  U = struct; % supplementary for construction 
  
  for f1 = this.ids % find dist metric crossover points 
    
    % subplot(numel(this.ids),1,find(f1 == this.ids))
    % cla, hold on, C = lines(max(7,nF));
    
    d_set = []; 
    d_idx = []; 
    for f2 = this.ids
      if f1 == f2, continue, end
 
        d = (splines.coeffs(1,:,f1)-splines.coeffs(1,:,f2)').^2 + ...
            (splines.coeffs(2,:,f1)-splines.coeffs(2,:,f2)').^2 ;
       [d,u] =  min(d);
        d_set = [d_set d']; %#ok<AGROW>
        d_idx = [d_idx u']; %#ok<AGROW>        
        % plot(d,'Color',C(f2,:))   
    end
    
    u_fid = this.ids(this.ids ~= f1);
    
   [d_min,d_fid] = min(d_set,[],2);    
    d_fid = u_fid(d_fid);    
    xop = find(d_fid ~= d_fid([2:end 1]))';    
    xop(d_min(xop) > 3*min(d_min(xop))) = []; 
    
    if ~any(xop)
      error('TODO one-boundary fascicle')
    end
    
    % plot(xlim, [2 2]*info.resol.^2,'k--')    
        
    U.dists{f1} = d_set;
    U.coeff{f1} = d_idx;
    U.xo_row{f1} = xop;    
    U.xo_dist{f1} = d_min(xop);    
  end
  
  clear d_fid d_idx d_min d_set f1 f2 sel u u_fid xod xop d
  
  %%
  U.xpi = {zeros(nF)}; 
  for f1 = this.ids % identify corner coefficient ids 
    %%
    if numel(U.xo_row{f1}) > 1
      error('TODO: untangle multiple corner points')
    end
    
    p_eff = zeros(size(U.dists));    
    for f2 = this.ids
      if f1 == f2, p_eff(f2) = U.xo_row{f1}; continue, end
      u_fid = this.ids(this.ids ~= f2);
      p_eff(f2) = U.coeff{f2}(U.xo_row{f2},u_fid == f1);
    end
    p_eff(p_eff == 0) = [];
        
    U.xpi{1}(f1,this.ids) = p_eff;    
        
    omo = min(U.dists{f1},[],2) > 3*min(U.xo_dist{f1});     
    omo = conv(omo([end 1:end 1]),[1;1;1],'valid') >= 2;
    
    u_fid = this.ids(this.ids ~= f1);

    for rr = find(omo ~= omo([2:end 1]))'
      
      if omo(rr), rr = rr-1; end %#ok<FXSET>
      % fprintf('f%d: rr=%d\n',f1,rr)
      [~,sel] = min(U.dists{f1}(rr,:)); 
      coe2 = U.coeff{f1}(rr,sel);
      
      mask = true(nF); % transpose mask
      mask([f1 u_fid(sel)],u_fid(sel)) = 0; 
      for pp = 2:numel(U.xpi)+1
        if pp > numel(U.xpi), U.xpi{pp} = zeros(nF);
        end
        if all(U.xpi{pp}(mask) == 0) && all(U.xpi{pp}(~mask) ~=0)
            break
        end
      end
      U.xpi{pp}(f1,f1) = rr;
      U.xpi{pp}(u_fid(sel),f1) = coe2;
    end
    
    
  end
  
  U.xco = zeros(numel(U.xpi),nF);  
  for pp = 1:numel(U.xpi) % collapse corners 
    xyz = []; 
    
    for f1 = this.ids
      if ~any(U.xpi{pp}(f1,:)), continue, end      
      coe = U.xpi{pp}(f1,U.xpi{pp}(f1,:)>0);
      U.xco(pp,f1) = median(coe);
      xyz = [xyz splines.coeffs(:,coe,f1)];
    end
    this.p_corner(pp,:) = mean(xyz,2);
  end
  
  for f1 = this.ids % order spline segments 
        
    inc = find(U.xco(:,f1)); 
    [~,ord] = sort(U.xco(inc,f1),'ascend');
    ord = inc(ord);
    this.s_order{f1} = []; 
    
    for rr = 1:numel(ord)
      r2 = mod(rr,numel(ord))+1;
      sel = all(this.s_corner == ord([rr r2])',2) | ...
            all(this.s_corner == ord([r2 rr])',2);
      if any(sel), 
        this.s_order{f1}(1,end+1) = find(sel,1); 
      else 
        this.s_corner = [this.s_corner; ord([rr r2])'];
        this.s_order{f1}(1,end+1) = size(this.s_corner,1);
      end
    end
  end

  U.object = {}; 
  
  U.xco(U.xco>0) = mod(U.xco(U.xco>0),size(splines.coeffs,2))+1; % tweak around
  
  for f1 = this.ids % get source spline object 
    nP = size(splines.coeffs,2); 
    U.object{f1} = spline(1:nP, splines.coeffs(:,:,f1)); 
  end
 
  %%
  if debug_plot, figure(33), clf, hold on, C = lines(7); end
  
  for pp = 1:size(this.s_corner,1) % compute segment interpolations 
    
    c1 = this.s_corner(pp,1);
    c2 = this.s_corner(pp,2);
    pair = (U.xco(c1,:).* U.xco(c2,:) ~= 0); 
    coe = U.xco([c1 c2],pair);     
    nP = max(3,round(mean(range(coe))));
    
    xyz = [];    
    for f1 = find(pair)
      
      % if (U.xco(c2,f1) > U.xco(c1,f1), f1 ~= find(pair,1))
        qp = linspace(U.xco(c1,f1),U.xco(c2,f1), nP);
        qp = conv(qp,[1 1]/2,'valid');
        
        if f1 == find(pair,1)          
             do_reverse = (sum(pair) == 1 && mean(diff(qp)) < 0); 
             if do_reverse
           px = size(splines.coeffs,2);
           qp = linspace(U.xco(c1,f1),U.xco(c2,f1)+px, nP);
           qp = conv(qp,[1 1]/2,'valid');
           qp = mod(qp-1,px+1);
             end
        elseif mean(diff(qp)) > 0, % do_reverse = true;
          
           px = size(splines.coeffs,2);
           qp = linspace(U.xco(c1,f1)+px,U.xco(c2,f1), nP);
           qp = conv(qp,[1 1]/2,'valid');
           qp = mod(qp-1,px+1);
           
          % do_reverse = mean(diff(qp)) > 1.2;
        end
      xyz = cat(3,xyz, ppval( U.object{f1}, qp));      
    end
    
    this.splines{pp,1} = mean(xyz,3); 
    
    xyf = [this.p_corner(this.s_corner(pp,1),:)'  this.splines{pp,1} ...
           this.p_corner(this.s_corner(pp,2),:)' ];

    del = ( xyf(:,2:end-1) - xyf(:, 1:end-2) );
    hat = ( xyf(:,3:end) - xyf(:, 1:end-2)   ); 
    hat = sum(hat .* del ) ./ sum( del.^2 ) ;
    
    if hat(1) < 1.25, this.splines{pp,1}(:,1) = []; end
    if hat(end) < 1.25, this.splines{pp,1}(:,end) = []; end    
    
    if debug_plot
      %%
      for q = 1:size(xyz,3)
        plot(xyz(1,:,q),xyz(2,:,q),'-','Color',C(q,:),'LineWidth',0.8,'userdata',qp)
        plot(xyz(1,1,q),xyz(2,1,q),'s','Color',C(q,:),'markersize',5,'userdata',qp)
        text(xyz(1,1,q)+0.005,xyz(2,1,q),num2str(pp),'Color',[.2 .2 .2],'FontSize',7)
      end
      plot(this.splines{pp,1}(1,:),this.splines{pp,1}(2,:),'Color',[.2 .2 .2],'LineWidth',1.2,'userdata',qp)
    end
  end
  
  if debug_plot
    plot(this.p_corner(:,1), this.p_corner(:,2), 'ko','MarkerFaceColor','k')
    axis image, tools.tidyPlot
  end
  
  for pp = 1:size(this.p_corner,1) % fix sharp corners in volume 
    
    [a,b] = find(this.s_corner == pp); 
    
    xyz = [cellfun(@(s) s(:,1:2), this.splines(a), 'unif', 0) ...
           cellfun(@(s) s(:,end-(0:1)), this.splines(a), 'unif', 0) ];
    xyz(b==2,1) =xyz(b==2,2); xyz(:,2) = [];
    xyz = cat(3,xyz{:}) - this.p_corner(pp,:)' ; 
    
    for ii = 1:size(xyz,3)
      if b(ii) == 1
           this.splines{a(ii)} = [ xyz(:,2,ii)*0.25 + this.p_corner(pp,:)' ...
                                   this.splines{a(ii)} ];
      else this.splines{a(ii)} = [ this.splines{a(ii)} xyz(:,2,ii)*0.25 + ...
                                   this.p_corner(pp,:)' ];
      end
    end
  end
  
  
  %%
  
  cs_id = numel(compound_structure);
  
  this.p_gmsh_ids = arrayfun(@(i) sprintf('p_c%d_%d',cs_id,i), ...
                                   1:size(this.p_corner,1),'unif',0)';
  
  this.s_gmsh_ids = arrayfun(@(i) sprintf('s_c%d_%d',cs_id,i), ...
                                   1:size(this.s_corner,1),'unif',0)';
  if isempty(compound_structure)
       compound_structure = this;
  else compound_structure(end+1) = this;
  end
  
  fprintf('Treating %s ids # %s%c as compound structure\n', name, ...
                    sprintf('%d,',this.ids), 8);

  clear cs_id rr r2 ans c1 c2 coe coe2 f1 f2 inc nF nP omo ord px 
  clear pp qp xyz mask sel pair q p_eff u_fid xyf del hat
end


%%
gmsh = ''; 
w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:});

if ~this.s_written
  %% write: corner points, boundary splines
  
  sstr = ''; 
  gmsh = w_(gmsh,'//Compound %s\np = newp;',name);
  for pp = 1:numel(this.p_gmsh_ids)    
    gmsh = w_(gmsh,'%s=p+%d;',this.p_gmsh_ids{pp},pp-1);    
  end
  
  for pp = 1:numel(this.p_gmsh_ids)    
    gmsh = w_(gmsh, 'Point(%s) = {%0.1f,%0.6f,%0.6f, lc} ;', ...
                          this.p_gmsh_ids{pp}, -splines.z, ...
                          this.p_corner(pp,2), this.p_corner(pp,1));
  end
    
  for ss = 1:numel(this.s_gmsh_ids)
  
    sstr = ''; 
    gmsh = w_(gmsh,'\n// for %s_%s\np = newp;',name,this.s_gmsh_ids{ss});

    for pp = 1:size(this.splines{ss},2) % last point is duplicate
        gmsh = w_(gmsh, 'Point(p+%d) = {%0.1f,%0.8f,%0.8f, lc} ;', ...
                                  pp-1, -splines.z, this.splines{ss}([2 1],pp));
        sstr = sprintf('%s,p+%d',sstr,pp-1);
    end
    sstr = strrep(sstr,',p+0','p');

    gmsh = w_(gmsh,'\nPrintf("Spline %s_%s : %%g - %%g", p, p+%d) ;', ...
                              name,this.s_gmsh_ids{ss},pp-1);
    gmsh = w_(gmsh,'%s = newreg;\nSpline(%s) = {%s,%s,%s} ;', ...
                this.s_gmsh_ids{[ss ss]}, ...   
                this.p_gmsh_ids{this.s_corner(ss,1)},  sstr, ...
                this.p_gmsh_ids{this.s_corner(ss,2)});
  end
  
  idx = cellfun(@(ids) isempty(setxor(this.ids,ids)), {compound_structure.ids});
  compound_structure(idx).s_written = true; 
end


sstr = sprintf(',%s',this.s_gmsh_ids{this.s_order{xy_id}});

gmsh = w_(gmsh,'\nPrintf("Loop %s_%02d : { %s }") ;', name, ... 
                                                xy_id, sstr(2:end));
gmsh = w_(gmsh,'e = newreg;\n Curve Loop(e) = {%s} ;', sstr(2:end));
gmsh = w_(gmsh,'Plane Surface(e+1) = e ;');
gmsh = w_(gmsh,'su[] = Extrude{%0.1f,0,0}{ Surface{e+1}; };', 2*splines.z);

return

%% Utilities
  
function is = in_loop(xy,loop) % are xy contained in the loop? 
% this is redundant with inPolygon I think

xy = xy';
loop = loop'; % Transpose to row-major

is = false(size(xy(:,1)));

for ii = 1:size(xy,1)
    
    dx = sign(loop(:,1) - xy(ii,1));
    cx = find(dx ~= circshift(dx,[1 0])); 
    if isempty(cx), continue, end
    cn = 0*cx; 
        
    % clf, hold on, C = lines(7);
    % plot(path(:,1),path(:,2),'color',[0 0 0 0.3])
    % plot(xy(ii,1),xy(ii,2),'s','LineWidth',1.2,'Color',C(2,:))
    % plot(path(cx,1),path(cx,2),'.k','MarkerSize',8)    
    % plot([1 1]*xy(ii,1),ylim,'--','LineWidth',1.2,'Color',C(2,:))
    
    for pp = reshape(cx,1,[])        
        % plot(path(pp-[0 1],1),path(pp-[0 1],2),'-k','LineWidth',1.2)
        
        if all(loop(pp-[0 1],2) >= xy(ii,2)), cn(cx==pp)=1; continue
        elseif all(loop(pp-[0 1],2) < xy(ii,2)),            continue
        end
        
        u = interp1(loop(pp-[0 1],1),loop(pp-[0 1],2),xy(ii,1));        
        cn(cx==pp) = (u > xy(ii,2));
    end
    
    is(ii) = mod(sum(cn),2)==1;
end
