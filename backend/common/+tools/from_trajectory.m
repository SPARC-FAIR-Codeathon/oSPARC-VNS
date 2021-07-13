

function [xyz,value] = from_trajectory(EM,nerve,xy, varargin)
% tools.from_trajectory( EM, fascicle, xy, [-f #], [-plot], [sensor] )
% 
% This code replicates what is done in mesh.insert_gmsh_fascicles. 
% IT's a tangled mess of things getting flipped this way and that
% becuase the XYZ coordinate system isn't super coherent between all the
% different software packages. 
% 
% the trajectory is assumed to be in mm, as are axon locations and
%   nerve.coeff. If nerve.coeff or axon xy are outside EM.info.DomainSize,
%   they are transformed from um to mm. 
% 
% You can also call this to determine if 3D trajectories are needed: 
%   is_3d = tools.from_trajectory( EM, fascicle )

if nargin == 2 
  xyz = isfield(nerve,'trajectory') || ( isfield(EM.info,'nerve') && ...
                     isfield(EM.info.nerve,'FascicleTrajectory')) || ...
                     isfield(EM.info,'FascicleTrajectory')        || ...
                     isfield(EM.info,'AxonTrajectory') ;
  return
end

basic_output = false; 

% Step 1: locate the relevent trajectories. 
if isfield(EM.info,'AxonTrajectory')
    error TODO_implement_axon_trajectory
elseif isfield(EM.info,'FascicleTrajectory'), 
    f_data = EM.info.FascicleTrajectory; 
elseif isfield(EM.info,'nerve') && isfield(EM.info.nerve,'FascicleTrajectory')
    f_data = EM.info.nerve.FascicleTrajectory;
elseif isfield(nerve,'trajectory'), f_data = nerve.trajectory; 
else basic_output = true; % no trajectory data
end

named = @(v) strncmpi(v,varargin,length(v));  
get_ = @(v) varargin{find(named(v))+1};

norm_ = @(v) v./sqrt(sum(v.^2,2)); % normal vector

ff = 1; 
if any(named('-f')), ff = get_('-f'); end

%%
coeff_XY = nerve.coeffs(:,:,ff)'; % Coeff ZYZ
coeff_XY = [0*coeff_XY(:,1) fliplr(coeff_XY)];
section_XY = [0*xy(:,1), fliplr(xy)]; % axon (or object) ZY intercept

if max(coeff_XY(:)) > max(EM.info.DomainSize), % convert um to mm
  coeff_XY = coeff_XY ./ 1e3; 
end
if max(section_XY(:)) > max(EM.info.DomainSize), 
  section_XY = section_XY ./ 1e3; 
end

if basic_output % kinda trivial

  f_XY = linspace(-EM.info.DomainSize(1),EM.info.DomainSize(1),101);
  f_XY = (f_XY' * [1 0 0] ) + mean(coeff_XY(2:end,[1 2 3]));
  
elseif iscell(f_data), f_XY = f_data{ff}; 
else                   f_XY = f_data;  
end

section_XY = section_XY - mean(coeff_XY(2:end,:)); % zero axon XY map
% section_XY = section_XY + f_XY(1,:); % translate to start point of trajectory

xyz = repmat(section_XY,1,1,size(f_XY,1));

for pp = 1:size(f_XY,1)-1

  rv = norm_(f_XY(pp+1,:)-f_XY(pp,:)); % rotate vector
  rot = [rv; norm_(cross(rv,[0 0 -1])); norm_(cross(rv,[0 1 0]))];
  xyz(:,:,pp) = section_XY * rot + f_XY(pp,:); 

end

xyz(:,:,pp+1) = section_XY * rot + f_XY(pp+1,:);

% xyz = xyz + f_XY(1,:) % translate to start point of trajector

do_VALUE = ~isempty(varargin) && isa(varargin{end},'scatteredInterpolant');

%%
if (nargout > 1 || any(named('-plot'))) && do_VALUE
  
  func = varargin{end}; 
  value = permute(func(xyz(:,3,:),xyz(:,2,:),xyz(:,1,:)),[1 3 2]);
 
else value = [];   
end





if any(named('-g')) % split into population groups 

    grp = get_('-g');
    nG = max(grp); 
    if any(named('-nG')), nG = get_('-nG'); end
    
    ax_xy = cell(nG,1);     
    for gg = 1:nG % For each group
      sel = grp == gg;
      ax_xy{gg,1} = xyz(sel,:,:);
    end
    if nargout == 1, xyz = ax_xy;
    else value = ax_xy;
    end
end




if ~any(named('-plot')), return, end
%%

cla reset

% 2D reference in XY plane
plot3(0*coeff_XY(1,:,ff),coeff_XY(2,:,ff),coeff_XY(1,:,ff),'k-')
hold on, axis equal, tools.tidyPlot
plot3(0*xy(:,1),xy(:,2),xy(:,1),'o')

fascicleName = sprintf('Fascicle%d',ff);
idx = EM.model.object_id{strcmpi(EM.model.object_name,fascicleName)};
% idx(EM.model.nodes(idx,1) < 0) = [];

plot3(f_XY(:,1),f_XY(:,2),f_XY(:,3),'-','LineWidth',1.2)

p.Vertices = EM.model.nodes(:,[1 2 3]);
p.Faces = [EM.model.elems(idx,1:3); EM.model.elems(idx,[1 3 4]);...
           EM.model.elems(idx,2:4); EM.model.elems(idx,[1 2 4])];

p.FaceColor = 'none';
p.EdgeColor = 'k';
p.EdgeAlpha = 0.2;
p.Clipping = 'off';
patch(p), tools.tidyPlot, grid on

C = lines(7);

for ii = 1:size(section_XY,1)

  pxyz = permute(xyz(ii,:,:),[3 2 1]);

  if do_VALUE
    patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),value(ii,:),'EdgeColor','interp','Marker','none','Clipping','off');
  else
    plot3(pxyz(:,1),pxyz(:,2),pxyz(:,3),'-','clipping','off','color',[C(2,:) 0.5])    
  end
end

axis auto, axis tight, colormap(tools.magma)
% axis([-0.1 0.1 -0.05 0.15 -0.2 0.2])
