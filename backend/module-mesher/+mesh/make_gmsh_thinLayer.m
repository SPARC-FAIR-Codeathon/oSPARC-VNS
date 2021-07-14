function em = make_gmsh_thinLayer(filename,varargin)
% Function em = make_gmsh_thinLayer(filename, ... )
%  reads a gmsh .geo file and add thin-layers, possibly overriding the GMSH
%  object names using "-name {thing:thing thing:thing}"

tools.setupEIDORS;

named = @(v) strncmpi(v,varargin,length(v)); 

%% Load mesh which needs thin-layers inserted 
if ~exist('filename','var'), 
    filename = tools.cache('PATH','pelvic_nerve.msh'); 
end
if strncmpi(fliplr(filename),'oeg.',4) % .geo extension
    filename = strrep(filename,'.geo','.msh'); 
end


if ~exist(filename,'file'), 
  error('could not find mesh file "%s"',filename)
end

fprintf('Reading "%s" ...\n', filename)

[~,nodes,~,object] = mesh.gmsh_read_mesh(filename,true); % <<<< EIDORS call

% WARNING to future users: as of May 2020, there's been a lot of work and 
% development of gmsh_read_mesh in EIDORS to handle the different file
% versions of the .msh file output by GMSH. I've got a hacked version of
% gmesh_read_mesh which I shared with the EIDORS community but there's no
% garuntee that future versions of GMSH or EIDORS will not crash here.
% 
% If you can't get this working, contact me to get my hacked
% gmsh_read_mesh, which worked with GMSH version 4.4.1 - CDE
%   $MeshFormat 4.1 0 8

while any(named('-n')) % Rename GMSH output objects if needed 
   disp('Overwriting physical volume names from gmsh file:')
   newname = varargin{find(named('-n'))+1};
   if any([newname{:}] == ':') % defines string-based renaming
       
       if any(cellfun(@(s) ~any(s==':'), newname)) || ...
                    numel(newname) < numel(object)
           cellfun(@disp, {object.name})
           error('all objects must be included in name list')
       end
       
       seq = cellfun(@(s) find(strncmpi({object.name}, s, ...
                                         find(s==':')-1)), newname);
       newname(seq) = regexprep(newname,'^[^:]+:','');
   end
   
   for ii = 1:numel(object)
     if ii > numel(newname) || isempty(newname{ii})
       fprintf('%2d: %8s -> [unchanged]\n',ii,object(ii).name)
       continue
     end
     fprintf('%2d: %8s -> %s\n',ii,object(ii).name, newname{ii})
     object(ii).name = newname{ii};
   end
   break
end

do_thin_layer = strncmpi({object.name},'P_',2);

if ~any(do_thin_layer)    
    warning('mk_thinLayer:nada',['No patch layers imported from gmsh ' ...
         '(looking for "P_[name]", where [name] is interior of patch)' ...
         ' try using -name to renamed these %d objects?'], numel(object))
end

if nargin > 1, layer_thickness = varargin{1}; 
else           layer_thickness = 1e-3; % 5e-3; 
end

%%
if 0 % Plot mesh to be thin-layered
  %%  
  clf, xyz = nodes.xyz; %#ok<UNRCH>
  trimesh(object(find(do_thin_layer,1)).nodes,xyz(:,1),xyz(:,2),xyz(:,3),'Clipping','off')
  axis equal off, hold on  
  if exist('ee','var')  
    axis(xyz(ee,[1 1 2 2 3 3]) +[-.1 .1 -.1 .1 -.1 .1])  
    plot3(xyz(ee,1),xyz(ee,2),xyz(ee,3),'rs','Clipping','off')
    plot3(xyz(eid,1),xyz(eid,2),xyz(eid,3),'r.','Clipping','off')  
    plot3(xyz(ee,1)+[0 0.02*dv(1)],xyz(ee,2)+[0 0.02*dv(2)],xyz(ee,3)+[0 0.02*dv(3)],'k-','Clipping','off')    
  end
  % if ei == nP, plot3(xyz0(:,1),xyz0(:,2),xyz0(:,3),'.','Clipping','off')
  % end
end

if all(layer_thickness == 0), do_thin_layer(:) = false; end
if numel(layer_thickness) == 1, 
  layer_thickness = layer_thickness * ones(1,sum(do_thin_layer)); 
end

% persistent DEBUG_PERM DEBUG_PID
% if isempty(DEBUG_PERM), DEBUG_PERM = flipud(perms(1:4)); DEBUG_PID = 1;
% else DEBUG_PID = DEBUG_PID + 1;
% end

compound = []; 

%%
for ff = find(do_thin_layer) % Generate thinlayers 
  %%
  
  if object(ff).dim == 3
      warning('mk_thinLayer:not2d', ...
              'Object %d (%s) has dimension %d (expected %d)' , ...
                      ff, object(ff).name, object(ff).dim,2)
      continue
  end
    
  if ~isempty(compound)
    error TODO_check_compound_strucutres
  end
    
  if isempty(object(ff).nodes) % GMSH doesn't always export correctly
    object = find_object_boundary(nodes, object, ff); 
  end
      
  if any(named('-c')) 
    [nodes, object, compound] = insert_compound_boundary(nodes,  ...
            object, compound, ff, do_thin_layer, layer_thickness); 
    error TODO_validate_multifascicle_thinlayer_mesh_sub-130
    
  else 
    % warning('a:b','DEV-MODE skipping normal fascicle %d', ff)
    % continue
    
    width = layer_thickness(find(do_thin_layer) == ff);    
    [nodes, object] = insert_simple_boundary(nodes, object, ff, width); 
  end

  if exist('DEBUG_PERM','var'), break,  end  
end

%% convert to EIDORS model struct and save 

clear ff node list enew xyz0 nP adj ei eid ee tt dv loc vol sel adj etet
em = construct_fwd_model(object,nodes); 


if exist('DEBUG_PERM','var')
    %%
    
    if ~exist('em','var')
        em = construct_fwd_model(object,nodes); 
    end
    
    clf,  show_fem(em)
    h = get(gca,'Children');
    h.EdgeAlpha = 0.5;
    h.FaceAlpha = 0.5;
    h.FaceColor = [1 1 1]; % [1 .8 .9];
    
    axis([4.8 5.01 0 0.2 0.05 0.25])
    
    set(gca,'CameraPosition',[ 6.3 -0.86 0.63])
    % axis([4.8 4.9 0.05 .2 -0 0.16])

    title(sprintf('[ %s]',sprintf('%d ',DEBUG_PERM(DEBUG_PID,:)))) %#ok<USENS>
    fprintf('[%cDEBUGGING THINLAYER %d/%d]%c\n', ...
                        8,DEBUG_PID,size(DEBUG_PERM,1),8)
    %%
    clear, return
end

filename = strrep(filename,'.msh','-thin.msh.mat'); 
fprintf('Saving "%s" ...\n', filename)
save(filename,'-struct','em')

return

%% CORE functions 
function object = find_object_boundary(nodes, object, index)


  %%
  inside = strcmp({object.name},object(index).name(3:end));    
  inside = object(inside);

  faces = inside.nodes(:,[1 2 3  1 2 4  1 3 4  2 3 4]);
  faces = [faces(:,1:3); faces(:,4:6); faces(:,7:9); faces(:,10:12)];
 [faces,~,fidx] = unique(sort(faces,2),'rows');    

  f_count = 0*faces(:,1);
  for ii = 1:size(fidx,1), f_count(fidx(ii)) = f_count(fidx(ii))+1; end

  faces(f_count > 1,:) = [];    
  x = reshape(nodes.xyz(faces,1),[],3);    
  xmax = max(nodes.xyz(:,1)); 
  
  domain_bound = all(abs(abs(x)-xmax) < 1e-9, 2); 
  faces(domain_bound,:) = []; 
    
%     clf
%     trimesh(faces,xyz(:,1),xyz(:,2),xyz(:,3),'FaceColor','none','EdgeColor','k')    
%     axis image, xlim([-1 1]*max(abs(ylim))-6)
    object(index).nodes = faces;

function [nodes, object] = insert_simple_boundary(nodes, object, index, width)

  ff = index; 
  nF = evalin('caller','sum(do_thin_layer)'); 

  xyz = nodes.xyz;
  node = object(ff).nodes';
  
  [list,~,enew] = unique(node);
  enew = reshape(enew,3,[]); 
  xyz0 = zeros(length(list),3);
  nP = length(list);
  
  %%
  
  adj = find(cellfun(@(e) any(ismember(e(:), object(ff).nodes(:))), ...
                                      {object.nodes}));
  adj(adj == ff) = []; % don't test this object
  interior = strcmp({object(adj).name},strrep(object(ff).name,'P_',''));
  
  if ~any(interior)
    warning('mk_thinLayer:noInterior', ...
            'Object %d (%s) does not have a defined inside (expected %s)' , ...
                    ff, object(ff).name, strrep(object(ff).name,'P_',''))
  end
  
  adj(interior) = [];  
    
  tic, tools.printInfo;   
  for ei = 1:nP % for each point in list 
      
    ee = list(ei);    
    if toc > 0.05, tic
      if isdeployed
          frac = 0.6*((ff-1)/nF + (ei/nP/nF));
          fprintf('progress: %0.3f%%\n', 100*frac + 30)
      else
        tools.printInfo('%s [%d/%d] (%0.2f%%)', object(ff).name, ei, nP, 100*ei/nP)
      end
    end
    
    eid = node(:,any(enew == ei,1));
    for tt = 1:size(eid,2) % circ-permute to align sel==ee
      if     eid(1,tt) == ee, eid(:,tt) = eid([2 3 1],tt);
      elseif eid(2,tt) == ee, eid(:,tt) = eid([3 1 2],tt);
      end
    end    
    loc = reshape(xyz(eid,:) - xyz(ee,:),3,[],3);
    % Make the cross-product normal for each triangle
    dv = [ loc(1,:,2).*loc(2,:,3) - loc(1,:,3).*loc(2,:,2); ...
           loc(1,:,3).*loc(2,:,1) - loc(1,:,1).*loc(2,:,3); ...
           loc(1,:,1).*loc(2,:,2) - loc(1,:,2).*loc(2,:,1)];
    dv = sum(dv,2)'; % take the 'average' and normalise 
    dv = dv / sqrt(sum(dv.^2));
    xyz0(ei,:) = xyz(ee,:) + width * dv; 
    
    for vol = adj % for each adjacent volume (excluding interior)
      object(vol).nodes(object(vol).nodes == ee) = ei + size(xyz,1); 
    end  
  end
  tools.printInfo('%s [%d/%d] (%0.2f%%)\n', object(ff).name, ei, nP, 100*ei/nP)

  enew = enew + size(xyz,1); % new node indices
  etet = zeros(4,3*size(enew,2)); % new simplexes  
  for tt = 1:size(enew,2) % build tet list from triangles    
    
     [~,idx] = sort(node(:,tt));
     
     enew_ins = [[node(:,tt); enew(idx(1),tt)] ...
                 [enew(:,tt); node(idx(3),tt)] ...
                 [node(idx(2:3),tt); enew(idx(1:2),tt)]];

    etet( 12*(tt-1)+1 : 12*tt ) = enew_ins;
  end
  
  % Insert new data into the mesh document structure
  object(ff).dim = 3;   
  object(ff).nodes = etet';
  
  sel = (nodes.obj(:,2) == object(ff).tag);
  nodes.xyz = [nodes.xyz; xyz0];     
  nodes.obj(sel,1) = 3; 

return


function [nodes, object, compound] = insert_compound_boundary(nodes,  ...
                 object, compound, index, do_thin_layer, layer_thickness)
          

this = object(index); 

this.basename = regexprep(this.name(3:end),'\d+','');
is_internal = strncmp({object.name},this.basename,length(this.basename));

searched = false(size(do_thin_layer)); 
shares_edge = false(size(do_thin_layer)); 
shares_edge(index) = 1; 

while any(shares_edge & ~searched) % search for touching objects 
  n = cat(1,object(shares_edge & ~searched).nodes);
  searched = searched | shares_edge;
  shares_edge = cellfun(@(x) any(ismember(n(:), x(:))), {object.nodes}); 
  shares_edge(~is_internal) = 0;
end

C.vol_idx = find(shares_edge);
C.surf_idx = cellfun(@(n) find(strcmp({object.name},['P_' n])), ...
                                     {object(C.vol_idx).name});
C.done = cell(size(C.vol_idx));

shares_edge = cellfun(@(x) any(ismember(n(:), x(:))), {object.nodes}); 
shares_edge([C.surf_idx C.vol_idx]) = 0;

for ii = C.surf_idx % get object boundaries (if needed)
  if isempty(object(ii).nodes)
    object = find_object_boundary(nodes, object, ii);
  end
  C.done{C.surf_idx == ii} = false(size(object(ii).nodes(:,1))); 
end

%% 

new_xyz = [];
new_tet = [];

tic, tools.printInfo;   
for ss = C.surf_idx
  for ff = 1:numel(C.done{ii}) % for each face in list    
    if C.done{ii}(ff), continue, end

    if toc > 0.05, tic
        tools.printInfo('%s [%d/%d] (%0.2f%%)', object(ff).name, ei, nP, 100*ei/nP)
    end

    
    error todo_facetype
    
  end
end


% for ei = 1:nP % for each point in list 
%       
%     ee = list(ei);
%     if toc > 0.05, tic
%         tools.printInfo('%s [%d/%d] (%0.2f%%)', object(ff).name, ei, nP, 100*ei/nP)
%     end
%     
%     eid = node(:,any(enew == ei,1));
%     for tt = 1:size(eid,2) % circ-permute to align sel==ee
%       if     eid(1,tt) == ee, eid(:,tt) = eid([2 3 1],tt);
%       elseif eid(2,tt) == ee, eid(:,tt) = eid([3 1 2],tt);
%       end
%     end    
%     loc = reshape(xyz(eid,:) - xyz(ee,:),3,[],3);
%     % Make the cross-product normal for each triangle
%     dv = [ loc(1,:,2).*loc(2,:,3) - loc(1,:,3).*loc(2,:,2); ...
%            loc(1,:,3).*loc(2,:,1) - loc(1,:,1).*loc(2,:,3); ...
%            loc(1,:,1).*loc(2,:,2) - loc(1,:,2).*loc(2,:,1)];
%     dv = sum(dv,2)'; % take the 'average' and normalise 
%     dv = dv / sqrt(sum(dv.^2));
%     xyz0(ei,:) = xyz(ee,:) + width * dv; 
%     
%     for vol = adj % for each adjacent volume (excluding interior)
%       object(vol).nodes(object(vol).nodes == ee) = ei + size(xyz,1); 
%     end  
%   end
%   tools.printInfo('%s [%d/%d] (%0.2f%%)\n', object(ff).name, ei, nP, 100*ei/nP)
% 
%   enew = enew + size(xyz,1); % new node indices
%   etet = zeros(4,3*size(enew,2)); % new simplexes  
%   for tt = 1:size(enew,2) % build tet list from triangles    
%     
%      [~,idx] = sort(node(:,tt));
%      
%      enew_ins = [[node(:,tt); enew(idx(1),tt)] ...
%                  [enew(:,tt); node(idx(3),tt)] ...
%                  [node(idx(2:3),tt); enew(idx(1:2),tt)]];
% 
%     etet( 12*(tt-1)+1 : 12*tt ) = enew_ins;
%   end
%   
%   % Insert new data into the mesh document structure
%   object(ff).dim = 3;   
%   object(ff).nodes = etet';
%   
%   % sel = (entity.obj_ids == object(ff).tag);
%   % entity.obj_dim(sel) = 4;
%   
%   nodes.xyz = [nodes.xyz; xyz0];     
%   nodes.obj(sel,1) = 3; 
% % nodes.vert; 

error TODO
% insert_compound_boundary


%% build fwd_model structure
function fwd_mdl = construct_fwd_model(object,nodes)

  object = object([object.dim] == 3);

  mdl.nodes    = nodes.xyz;
  mdl.elems    = uint32(cat(1,object.nodes));
  mdl.boundary = [];
  mdl.boundary_numbers = []; 
  mdl.gnd_node = 1;
  mdl.name = 'Fascicle simulation';

  % Model Stimulation - nil for now. 
  % if ~isempty(stim_pattern)
  %     mdl.stimulation= stim_pattern;
  % end

  % Electrodes
  sel = find(arrayfun(@(x)strncmp(x.name,'Elec',4),object));    
  for ee = 1:length(sel), id = sel(ee);         
      this.name = object(id).name;
      this.nodes = unique(object(id).nodes); 
      this.z_contact = 1e-12; 
      mdl.electrode(ee) = this;
  end

  mdl.solve =      'eidors_default';
  mdl.jacobian =   'eidors_default';
  mdl.system_mat = 'eidors_default';

  fwd_mdl = eidors_obj('fwd_model', mdl);

 [index, order] = mk_mat_indices( object );
  fwd_mdl.boundary = find_boundary(fwd_mdl);
  fwd_mdl.object_name = {object(order).name};
  fwd_mdl.object_id = index;

return
    
%% Output cell array of indices into each material type
%    array order is sorted by length of material type
function [mat_indices, order] = mk_mat_indices( object )

    sel = find([object.dim] == 3);
    
    if any(cellfun(@isempty,{object(sel).nodes}))
        bad = cellfun(@isempty,{object.nodes});
        bad = sprintf(',%s',object(bad).names);
        
        error('Empty physical domains: {%s}. Usually this indicates %s', ...
               bad(2:end),'a serious issue with the model geometry.')
        
    end
    
    % & ~cellfun(@isempty,{object.nodes}));  ... if this is true there may
    % be other bad issues with the gmsh code? 
    
    index = arrayfun(@(e) 0*object(e).nodes(:,1)+e, sel,'Unif',0);
    index = cat(1,index{:});  

    n_e = cellfun(@length,{object(sel).nodes}); 
    [~,order] = sort(n_e,'descend'); % reverse sort

    mat_indices = cell(1, length(order));
    for ii = 1:length(order)
        mat_indices{ii}= find(index == order(ii));
    end