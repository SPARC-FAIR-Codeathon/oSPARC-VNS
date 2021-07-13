function[srf,vtx,fc,bc,simp,edg,mat_ind,obj_names] = gmsh_read_mesh(filename,raw_out)
%[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(filename)
% Function to read in a mesh model from Gmsh and saves it in
% five arrays; surface (srf), veritices (vtx), face no. (fc)
% volume (simp) and edges (edg)
%
% srf        = The surfaces indices into vtx
% simp       = The volume indices into vtx
% vtx        = The vertices matrix
% fc         = A one column matrix containing the face numbers
% edg        = Edge segment information
% filename   = Name of file containing NetGen .vol information
% mat_ind    = Material index
% phys_names = Structure of "Physical" entities in the mesh
%              .dim   = dimension
%              .name  = name (string)
%              .tag   = physical tag
%              .nodes = N-x-dim array of indices into vtx

% $Id: gmsh_read_mesh.m 3260 2012-06-30 14:40:10Z bgrychtol $
% (C) 2009 Bartosz Sawicki. Licensed under GPL V2
% Modified by Calvin Eiber <ceiber@unimelb.edu.au> 
%             James Snyder <jbsnyder@fanplastic.org>

if nargin == 0, % Because things should work nicely when you press "run"
  filename = dir('*.msh'); 
  if numel(filename) == 1, filename = filename.name; 
  else                     [fn,fp] = uigetfile('*.msh'); 
                           filename = [fp fn];
  end
end

fid = fopen(filename,'r');
obj_names = [];
entities = []; 

while 1
  tline = fgetl(fid);
  if ~ischar(tline); fclose(fid); break; end
  
  switch(tline)
    case '$Elements',      elements = parse_elements( fid );
    case '$Nodes',         nodes    = parse_msh_nodes( fid );
    case '$Entities',      entities = parse_geo_entities( fid );
    case '$PhysicalNames', obj_names = parse_names( fid );
  end
end

if ~isempty(obj_names)
  for ii = 1:numel(obj_names)
    % entities defines a bijective relationship between mesh entites and
    % physical groups. 
    these = cellfun(@(v) any(ismember(v,obj_names(ii).tag)), ...
                                        entities.phy_ids);
    if ~any(these), continue, end
        
    obj_names(ii).tag = entities.obj_ids(these); %#ok<AGROW>
    these = ismember([elements.geom_tag], obj_names(ii).tag);  
    obj_names(ii).nodes = cat(1,elements(these).simp); %#ok<AGROW>
  end
end

edg = [];
bc = [];

% Select 2d vs 3d model by checking if Z is all the same
if length( unique( nodes.xyz(:,3) ) ) > 1 
    vtx = nodes.xyz;
    % Type 2: 3-node triangle
    tri = ([elements.type] == 2);
    srf = cat(1,elements(tri).simp);
    % Type 3: 4-node tetrahedron
    tet = ([elements.type] == 3);
    simp = cat(1,elements(tet).simp);
    
    mat_ind = arrayfun(@(e) 0*elements(e).simp(:,1)+e, find(tet),'Unif',0);
    mat_ind = cat(1,mat_ind{:});
else
    vtx = nodes(:,2:3);
    tri = ([elements.type] == 2);
    simp = cat(1,elements(tri).simp);
    srf = [];
        
    mat_ind = arrayfun(@(e) 0*elements(e).simp(:,1)+e, find(tri),'Unif',0);
    mat_ind = cat(1,mat_ind{:});
end

elemtags = cat(1,elements.geom_tag);
fc = elemtags(tri,1);

if nargin < 2 || ~raw_out, return, end

srf = elements; 
vtx = nodes; 
fc = entities; 
bc = obj_names; 



    
    

end




function geo = parse_geo_entities( fid ) % $Entities

tline = fgetl(fid);
n_obj = sscanf(tline,'%d'); % nPoints, nCurves, nSurfaces, nVolumes

geo.obj_ids = [];
geo.obj_dim = []; 
geo.phy_ids = [];

for ty = 1:4 
  
  pids = zeros( n_obj(ty), 1 );
  tags = cell( n_obj(ty), 1 );

  for ii = 1:n_obj(ty)
  
    tline = sscanf(fgetl(fid),'%f');
    
    if ty == 1, tline(2:4) = []; 
    else        tline(2:7) = []; 
    end
    
    pids(ii) = tline(1); 
    tags{ii} = tline(2 + (1:tline(2)));
    
  end
  
  geo.obj_ids = [geo.obj_ids; pids]; 
  geo.obj_dim = [geo.obj_dim; 0*pids + ty];
  geo.phy_ids = [geo.phy_ids; tags]; 
  
end

end

function node = parse_msh_nodes( fid )
% Line Format: - this is OLD, CE 31/07/2019 
% node-number x-coord y-coord z-coord

% BEGIN_PTR = ftell(fid); 
% fseek(fid,BEGIN_PTR,-1)

n_obj = fscanf(fid,'%d',[1 4]);

node.obj  = zeros(n_obj(1),4);
node.vert = cell(n_obj(1),1);
node.xyz = zeros(n_obj(2),3);

% raw = fscanf(fid,'%f',[4 (n_obj(1)+n_obj(2))]); 
% ptr = 0; 

for id = 1:n_obj(1)

    node.obj(id,:) = fscanf(fid,'%d',[1 4]);
    n = node.obj(id,4); 

    node.vert{id} = fscanf(fid,'%d',[1 n]);
    xyz = fscanf(fid,'%f %f %f',[3 n]);
    
    node.xyz(node.vert{id},:) = xyz'; 
end

% mat = [[node.vert{:}]' node.xyz];

end

function names = parse_names( fid )
% Line Format:
% physical-dimension physical-number "physical-name"
tline = fgetl(fid);
n_rows = sscanf(tline,'%d');
names = struct('tag',{},'dim',{},'name',{});
for i = 1:n_rows
    tline = fgetl(fid);
    if exist('OCTAVE_VERSION'),  %#ok<EXIST>
         parts = strsplit(tline,' ');
    else parts = regexp(tline,' ','split');
    end
    nsz = size(names,2)+1;
    names(nsz).dim = str2double( parts(1) );
    names(nsz).tag = str2double( parts(2) );
    tname = parts(3);
    names(nsz).name = strrep(tname{1},'"','');
end
end

function elements = parse_elements( fid )
% Line Format:
% num-entities num-elements _?=1_ _?=num-elements_ 

n_rows = fscanf(fid,'%d',[1 4]);
n_rows = n_rows(1); % We'll grow the elements arrays as needed

% elements = struct('simp',{},'phys_tag',{},'geom_tag',{});
elements(n_rows).simp = [];
elements(n_rows).phys_tag = [];
elements(n_rows).geom_tag = [];
elements(n_rows).type = [];

for i = 1:n_rows(1)
    
    n = fscanf(fid,'%d',[1 4]);
    % Line Format: <tags = n-dims  obj-id  parent-id  n-simplices>
    
    % By default, first tag is the type of mesh (1D = 0, 2D = 1, ...) 
    % second is the number of parent elementary geometrical entity
    % third is parent physical entity
    % fourth is number of parent mesh partitions followed by
    % partition vertex list
    
    elements(i).type = n(1); 
    elements(i).phys_tag = n(3); 
    elements(i).geom_tag = n(2);     
    elements(i).simp = fscanf(fid,'%d',[2+n(1) n(4)])';
    elements(i).enum = elements(i).simp(:,1);
    elements(i).simp(:,1) = []; % First # in each row is elem number
end
end
