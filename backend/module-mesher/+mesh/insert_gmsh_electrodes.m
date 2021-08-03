
function gmsh = insert_gmsh_electrodes(varargin)

% named = @(v) strncmpi(v,varargin,length(v)); 

if nargin == 0, G = get_layout([]); 
elseif nargin >= 1 && isstruct(varargin{1})
     G = get_layout(varargin{1},varargin{2:end});
elseif nargin >= 1 && ischar(varargin{1})
     G = get_layout(load_layout(varargin{1}), varargin{2:end});
else G = get_layout([],varargin{:});
end

nE = size(G.EP,1);

if numel(G.ETI) < nE, G.ETI = repmat(G.ETI(:),nE,1); end
if numel(G.ID) < max(G.ETI), G.ID = repmat(G.ID,nE,1); end

% Assuming: Volume(1) is the electrode carrier (gmsh_CarrierIndex)
w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:}); 
gmsh = sprintf('\n// Electrodes\nel = newv;');

% Electrode cut-out in carrier
xyz = @(e) [ G.EP(e,:)-G.ED(G.ETI(e),:)./[2 1 2]-[0 G.ID(G.ETI(e)) 0] ...
                       G.ED(G.ETI(e),:)+[0 G.ID(G.ETI(e)) 0] ];

                   
for ee = 1:nE % Loop 1: electrode cut-out
    
    if strncmpi(G.EK{G.ETI(ee)}, 'circum',6)
        
      r = G.ID(G.ETI(ee))+G.ED(G.ETI(ee),2);
      if isfield(G,'cuff_ID'), r = r+G.cuff_ID; 
      elseif isfield(G,'carrier') && isfield(G.carrier,'cuff_IDr')
           r = r + G.carrier.cuff_IDr;      
      end
      if size(G.ED,2) == 2, a = 2*pi; 
      else a = G.ED(G.ETI(ee),3)/(r-G.ED(G.ETI(ee),2));
      end, a = min(a,2*pi); 
          
      xra = [G.EP(ee,1)-G.ED(G.ETI(ee),1)/2   G.ED(G.ETI(ee),1) r a];
      gmsh = w_(gmsh,'Cylinder(el) = {%0.3f,0,0,%0.3f,0,0,%0.3f,%0.9f};',xra);      
      gmsh = w_(gmsh,'Rotate{ {1,0,0}, {0,0,0}, %0.9f } { Volume{el}; }', pi/2-a/2);

    elseif ~isempty(G.ER) && G.ER(G.ETI(ee)) > 0 || ...
            strncmpi(G.EK{G.ETI(ee)}, 'ci',2)
      z = G.ID(G.ETI(ee))+G.ED(G.ETI(ee),2); 
      xyz_r = [G.EP(ee,:) - [0 z 0] [0 z 0] G.ER(G.ETI(ee))];
      gmsh = w_(gmsh,'Cylinder(el) = {%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,2*Pi};',xyz_r);
    else
      gmsh = w_(gmsh,'Box(el) = {%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f};', xyz(ee));
    end
    
    if ~isempty(G.EA) % ElectrodeAngle        
      gmsh = w_(gmsh,'Rotate{ {1,0,0}, {0,0,0}, %0.9f } { Volume{el}; }', deg2rad(G.EA(ee)));
    end
    gmsh = w_(gmsh,'BooleanDifference{ Volume{%d}; Delete; }{ Volume{el}; Delete; }', ...
                    G.gmsh_CI);
end


gmsh = w_(gmsh,'\n');

% Electrode volume
xyz = @(e) [ G.EP(e,:)-G.ED(G.ETI(e),:)./[2 1 2]-[0 G.ID(G.ETI(e)) 0] ...
                       G.ED(G.ETI(e),:) ];
                   
for ee = 1:nE % Loop 2: electrode objects
    
    if strncmpi(G.EK{G.ETI(ee)}, 'circum',6)
        
      r = G.ID(G.ETI(ee))+G.ED(G.ETI(ee),2);
      if isfield(G,'cuff_ID'), r = r+G.cuff_ID; 
      elseif isfield(G,'carrier') && isfield(G.carrier,'cuff_IDr')
           r = r + G.carrier.cuff_IDr;      
      end
      if size(G.ED,2) == 2, a = 2*pi; 
      else a = G.ED(G.ETI(ee),3)/(r-G.ED(G.ETI(ee),2));
      end, a = min(a,2*pi); 
          
      x = [G.EP(ee,1)-G.ED(G.ETI(ee),1)/2   G.ED(G.ETI(ee),1)];
      gmsh = w_(gmsh,'Cylinder(el+%d) = {%0.3f,0,0,%0.3f,0,0,%0.3f,%0.9f};',ee-1,x,r,a);      
      gmsh = w_(gmsh,'Rotate{ {1,0,0}, {0,0,0}, %0.9f } { Volume{el+%d}; }', pi/2-a/2, ee-1);
      
      r = r - G.ED(G.ETI(ee),2);
      gmsh = w_(gmsh,'Cylinder(el+%d) = {%0.3f,0,0,%0.3f,0,0,%0.3f,%0.9f};',ee,x,r,a);      
      gmsh = w_(gmsh,'Rotate{ {1,0,0}, {0,0,0}, %0.9f } { Volume{el+%d}; }', pi/2-a/2, ee);
      gmsh = w_(gmsh,'BooleanDifference{ Volume{el+%d}; Delete; }{ Volume{el+%d}; Delete; }', ...
                    ee-[1 0]);
                
    elseif ~isempty(G.ER) && G.ER(G.ETI(ee)) > 0 || ...
            strncmpi(G.EK{G.ETI(ee)}, 'ci',2)
      z = G.ID(G.ETI(ee))+G.ED(G.ETI(ee),2); 
      xyz_r = [G.EP(ee,:) - [0 z 0] [0 z-G.ID(G.ETI(ee)) 0] G.ER(G.ETI(ee))];
      gmsh = w_(gmsh,'Cylinder(el+%d) = {%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,2*Pi};',ee-1,xyz_r);
    else    
      gmsh = w_(gmsh,'Box(el+%d) = {%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f};', ee-1, xyz(ee));
    end
    
    if ~isempty(G.EA) % ElectrodeAngle
      gmsh = w_(gmsh,'Rotate{ {1,0,0}, {0,0,0}, %0.9f } { Volume{el+%d}; }', deg2rad(G.EA(ee)), ee-1);
    end
end

% if G.ER, G.ESR, or G.IR were non-zero we'd have to handle calling
% "Fillet{}{}{}" here as well as finding the relevent indices
% asset(isempty(G.ER),'TODO: implement electrode radius')
% assert(isempty(G.ESR),'TODO: implement electrode surface radius')
% assert(isempty(G.IR),'TODO: implement electrode inset radius')


gmsh = insert_domain_medium(gmsh,G,nE);


for ee = 1:nE % Loop 3: BooleanDifference
    gmsh = w_(gmsh,'Physical Volume("Elec%d") = {el+%d};', ee, ee-1);
end

return



function gmsh = insert_domain_medium(gmsh, G, nE)

w_ = @(g,f,varargin) sprintf(['%s\n' f],g,varargin{:}); 

layers = []; 

if any(G.DS < 0)
  nb = G.DS(G.DS < 0); 
  G.DS(G.DS < 0) = []; 
 switch sum(G.DS < 0)
  case 1, OB = [-G.DS(1) G.DS(1) -G.DS(2) G.DS(2)    nb    G.DS(3)];
  case 2, OB = [   nb(1) G.DS(1)    nb(2) G.DS(2) -G.DS(3) G.DS(3)];
  case 3, OB = [   nb(1) G.DS(1)    nb(2) G.DS(2)    nb(3) G.DS(3)];
  otherwise
     error('unable to handle %d negative bounds on G.DomainSize', sum(G.DS < 0))
 end 
else      OB = [-G.DS(1) G.DS(1) -G.DS(2) G.DS(2) -G.DS(3) G.DS(3)];
end
layers = G.DS(4:end);

if isfield(G,'DN') && ismember(numel(G.DN), numel(layers)+[0 1 2])
  layer_names = G.DN;
elseif isfield(G,'D') && ismember(numel(G.D), numel(layers)+[0 1 2])
  layer_names = G.D;
end

if numel(layer_names) <= numel(layers)
  for ii = (numel(layer_names)+1) : numel(layers)
    layer_names{ii} = sprintf('Tissue_%02d',ii); 
  end
end

layers = [layers OB(6)];

% if ~any(named('-interstitial')), return, end
gmsh = w_(gmsh,'\nvol_id = newv;'); % Add domain
for ii = 1:numel(layers)

  if ii == 1, z = [OB(5) layers(ii)-OB(5)];
  else z = [layers(ii-1) layers(ii)-layers(ii-1)];
  end
    
  gmsh = w_(gmsh,'Box(vol_id+%d) = {%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f}; // %s', ...
                             ii-1,  OB(1), z(1),OB(3), ...
                                    OB(2)-OB(1),z(2),OB(4)-OB(3), ...
                                                  layer_names{ii});
end

gmsh = w_(gmsh,'BooleanDifference{ Volume{vol_id}; Delete; }{ Volume{%d}; }', G.gmsh_CI);

for ee = 1:nE % Loop 3: BooleanDifference
    gmsh = w_(gmsh,'BooleanDifference{ Volume{vol_id}; Delete; }{ Volume{el+%d}; }', ee-1);
end

gmsh = w_(gmsh,'\nPhysical Volume("PDMS") = {%d};', G.gmsh_CI);
for ii = 1:numel(layers)
    gmsh = w_(gmsh,'Physical Volume("%s") = {vol_id+%d};',layer_names{ii},ii-1);
end


% parse varargin for geometry properties
function geom = get_layout(geom,varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if isempty(geom), geom = struct; end

if isfield(geom,'array') % extract 'array' and 'mesh' as 'geom'
  if isfield(geom,'mesh')
    g = geom.mesh;
    for f = fieldnames(geom.array)'
      g.(f{1}) = geom.array.(f{1});
    end
    geom = g; 
  else geom = geom.array; 
  end
end   

% Set up parameter defaults
default.ElectrodePositions = [1.85 0 0; 1.1 0 0; -1.1 0 0; -1.85 0 0];
default.ElectrodeDimensions = [0.2 0.1 1.8]; % L W H
default.ElectrodeRadius = []; 
% default.ElectrodeSurfaceRadius = [];
default.ElectrodeAngle = []; 
default.ElectrodeKind = {'Rectangular'};
default.ElectrodeTypeIndex = 1; 
default.InsetDepth = 0.1;
default.InsetRadius = 0;

default.DomainSize = [6 6 6]; % mm
default.DomainName = {'Interstitial'};

default.gmsh_CarrierIndex = 1; 

if any(named('-import'))
  if exist(get_('-import'),'file')
       geom = load(get_('-import'));
  elseif exist(tools.cache('path',get_('-import')),'file')    
       geom = load(tools.cache('path',get_('-import')));
  else warning('insert_gmsh:filenotfound', ...
               '%s "%s" not found. Using defaults...', ...
               'Electrode definition file', get_('-import'))
  end
end

for f = fieldnames(default)' % For each field    
  if ~isfield(geom,f{1}), geom.(f{1}) = default.(f{1}); end
  if any(named(f{1})), geom.(f{1}) = get_(f{1}); end
  if ~strncmpi(f,'Electrode',9), continue, end
  f_sh = strrep(f{1},'Electrode',''); % "Type" refers to "ElectrodeType", etc.
  if any(named(f_sh)), geom.(f{1}) = get_(f_sh); end    
end

%% Run through verification checks
geom = compress_fieldnames(geom); % removes lowercase letters

if ~iscell(geom.EK), geom.EK = {geom.EK}; end    
nE = size(geom.EP,1);
nT = max(geom.ETI); 

if numel(geom.ETI) < nE, geom.ETI = repmat(geom.ETI(:),nE,1); end
if numel(geom.EK) < nT,  geom.EK = repmat(geom.EK,nT,1); end
if numel(geom.ID) < nT,  geom.ID = repmat(geom.ID,nT,1); end

return


function geom = load_layout(filename)

has_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if any( filename == '~'), filename = tools.file(filename); end
      
if ~exist(filename,'file') && exist(tools.cache('path',filename),'file')
    filename = tools.cache('path',filename); 
end

if has_ext_(filename,'.mat') % from mat file
    % printf('Loading %s', source_file);
    geom = load(filename);  
    geom.FILENAME = filename;
elseif has_ext_(filename,'.xml') % from XML
    % printf('Loading %s', source_file);
    geom = tools.parseXML(filename);
    geom.FILENAME = filename;
elseif has_ext_(filename,'.json')    
    % printf('Loading %s', source_file);
    geom = tools.parse_json(filename);
    geom.FILENAME = filename;

else error('unknown filetype on file "%s", expected {.mat, .xml, .json}', filename)
end

function S = compress_fieldnames(S) 
% compress_fieldnames is a convenience function to turn verbose field names 
%  (good for saving data and as input/output arguments) into compact field
%  names (good for working with data, particularly geometry) by removing
%  lowercase letters. 

F_long = fieldnames(S);
F_short = regexprep(F_long,'[a-z]','');
for ii = 1:length(F_long)    
  if any(strcmp(F_short{ii},F_short(1:ii-1))), continue, end
  if isempty(F_short{ii}), continue, end
  if F_short{ii}(1) == '_', F_short{ii} = ...
      [regexp(F_long{ii},'^[^_]+_','match','once') F_short{ii}(2:end)];
  end
  S = renameStructField(S,F_long{ii},F_short{ii});
  
  % compress {} to [], clean up tools.parse_json mess
  if iscell(S.(F_short{ii})) && all(cellfun(@numel,S.(F_short{ii})) == 1)
      S.(F_short{ii}) = [S.(F_short{ii}){:}]; 
  end
end

if isfield(S,'EP') && iscell(S.EP)
    S.EP = cellfun(@(x) [x{:}], S.EP,'unif',0);
    S.EP = cat(1,S.EP{:});
end

function str = renameStructField(str, oldFieldName, newFieldName)
%RENAMESTRUCTFIELD renames oldFieldName to newFieldName in struct str
%
%   STR = RENAMESTRUCTFIELD(STR, OLDFIELDNAME, NEWFIELDNAME)
%   STR is the struct in which to rename the field
%   OLDFIELDNAME is the name of the field to rename
%   NEWFIELDNAME is the name to rename OLDFIELDNAME to
%

%   Copyright 2013-2014 The MathWorks, Inc.
if ~strcmp(oldFieldName, newFieldName)
    allNames = fieldnames(str);
    % Is the user renaming one field to be the name of another field?
    % Remember this.
    isOverwriting = ~isempty(find(strcmp(allNames, newFieldName), 1));
    matchingIndex = find(strcmp(allNames, oldFieldName));
    if ~isempty(matchingIndex)
        allNames{matchingIndex(1)} = newFieldName;
        str.(newFieldName) = str.(oldFieldName);
        str = rmfield(str, oldFieldName);
        if (~isOverwriting)
            % Do not attempt to reorder if we've reduced the number
            % of fields.  Bad things will result.  Let it go.
            str = orderfields(str, allNames);
        end
    end
end

