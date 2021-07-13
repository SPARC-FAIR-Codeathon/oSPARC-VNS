

function q = meshQuality(file)
% Compute mesh quality metrics for EIDORS mesh
% 
% Element quality was assessed according to the metrics proposed by 
% Joe and Liu \cite{Liu1994} with minimum element quality~0.527 and 
% average quality~0.887

if nargin == 0 || ~isfile(file)
  file = tools.file('~\source\mesh\pelvic_nerve-thin.msh.mat');
end


tools.setupEIDORS; 

m = load(file);
if isfield(m,'model'), m = m.model; end

q = calc_mesh_quality(m);
[~,m.name] = fileparts(file);

%% And display

fprintf('\n\n%s: %d elements\n==== MESH QUALITY METRICS ====\n', m.name, size(m.elems,1))

ok = cat(1,m.object_id{~strncmp(m.object_name,'P_Fascicle',3)});
pf = cat(1,m.object_id{ strncmp(m.object_name,'P_Fascicle',3)});

for ff = fieldnames(q.tet)'    
    val = q.tet.(ff{1})(ok,:);
    val = [min(val) mean(val)];
    fprintf('      min %s %0.3f, mean %0.3f\n',ff{1},val)
    q.summary.(['obj_' ff{1}]) = val;
    
    val = q.tet.(ff{1})(pf,:); 
    val = [min(val) mean(val)];   
    fprintf('ThinLayer %s %0.3f, mean %0.3f\n',ff{1},val)
    q.summary.(['P_' ff{1}]) = val;
end

q.name = m.name; 
q.n_elements = size(m.elems,1); 
q.n_verts = size(m.nodes,1); 

if nargout == 0, clear, end


