

function view_mesh(filename,varargin)
% function plots.view_mesh(filename,...)
% -prompt uses a file select dialog (default: read from ~/source/mesh/)
% -fascicle shows fascicles in blue
% -animate animates each named object from GMSH
% -no-elec hides electrodes (if relevent)

if nargin > 0 && ischar(filename) && ~isempty(filename) && filename(1) == '-'
  varargin = [{filename} varargin]; filename = ''; 
end

named = @(v) strncmpi(v,varargin,length(v)); 
tools.setupEIDORS; 

if nargin == 0 || isempty(filename) % output file
  filename = tools.file('out~/*.mat','-prompt');
elseif isstruct(filename), m = filename; filename = ''; 
elseif strcmpi(filename,'SOURCE')
    filename = tools.file('array~/pelvic_nerve.geo');
end

if exist('m','var') % arg-in supplied
elseif strncmpi(fliplr(filename),'oeg.',4) % .GEO file

  m = load(regexprep(filename,'\..*$','-thin.msh.mat')); 
  m.boundary_numbers = ones(size(m.boundary(:,1)));
  
elseif strncmpi(fliplr(filename),'tam.',4) % .MAT file
 
  m = load(filename);
  if isfield(m,'model'), m = m.model; end
  if ~isfield(m,'boundary_numbers') 
    m.boundary_numbers = ones(size(m.boundary(:,1)));
  end
  
else error('%s: unknown filetype for plots.%s',filename,mfilename)
end

if any(named('-no-e')) || any(named('-e')),  m.electrode = []; end

%% Show the volume and electrode pattern

cla reset, show_fem(m)
sel = cat(1,m.object_id{strncmp(m.object_name,'Elec',4)});
patch('Faces',m.elems(sel,:),'Vertices',m.nodes,'EdgeColor',[1 .4 .2], ...
                        'FaceColor','w','FaceAlpha',0.2,'EdgeAlpha',0.5,'userdata','elec')
h = get(gca,'Children');
h(end).EdgeAlpha = 0.2;
pause(0.05), hold on

tet2tri = @(n) [n(:,1:3); n(:,[1 2 4]); n(:,[1 3 4]); n(:,2:4)];


if any(named('-f'))

  blue = [.3 .5 .9];
  sel = cat(1,m.object_id{strncmp(m.object_name,'P_Fascic',4)});
  if isempty(sel)
    sel = cat(1,m.object_id{strncmp(m.object_name,'Fascic',4)});
    patch('Faces',tet2tri(m.elems(sel,:)),'Vertices',m.nodes,'EdgeColor',blue, ...
                            'FaceColor',blue,'FaceAlpha',0.4,'EdgeAlpha',0.4,'userdata','fasc')
  else
    patch('Faces',tet2tri(m.elems(sel,:)),'Vertices',m.nodes,'EdgeColor',blue, ...
                            'FaceColor',blue,'FaceAlpha',0.4,'EdgeAlpha',0.4,'userdata','fasc')
  end
  
end

if ~(any(named('-no-p')) || any(named('-nop')))

  sel = cat(1,m.object_id{strncmp(m.object_name,'PDMS',4)});
  patch('Faces',tet2tri(m.elems(sel,:)),'Vertices',m.nodes,'EdgeColor',[1 .96 .7], ...
                          'FaceColor','w','FaceAlpha',0,'EdgeAlpha',0.4,'userdata','pdms')
end

axis(axis)
% plot3([-1 0 0 0 0]+5,[0 0 1 0 0]-2,[0 0 0 0 1]-6.5,'Color',[.3 .3 .3],'LineWidth',1.1,'clipping','off')

set(gca,'CameraPosition',get(gca,'CameraPosition') .* [-1 -1 1])
hold on, axis off, 
title(strrep(filename,tools.file,'~'),'interpreter','none')



if any(named('-xy')),     arrange_axes_xy; 
elseif any(named('-yz')), arrange_axes_yz; 
end



%%

if ~any(named('-a')), return, end

for ii = 3:numel(m.object_id)     
    % be careful - there's two sets of indices for different objects which
    % aren't necessarially in the right order. 
    oid = m.object_id{ii}; %%#ok<UNRCH>
    h(1).Faces = tet2tri(m.elems(oid,:)); 
    title(m.object_name{ii})
    [~,~,b] = ginput(1); 
    if b == 27, return, end
end

clear oid i h b t



function arrange_axes_yz

  view([1 0 0])
  axis([-6.1 -5.95 0 0.3 -0.3 0.3])
  axis on

  h = get(gca,'Children');

  h(2).FaceAlpha = 0.1;
  h(2).LineWidth = 1.2;
  h(2).EdgeAlpha = 1;
  x = h(2).Vertices(:,1);

  keep = all ( x(h(2).Faces) < -6 + 10^-4 , 2); 
  h(2).Faces = h(2).Faces(keep,:);

  keep = all ( x(h(1).Faces) < -6 + 10^-4 , 2); 
  h(1).Faces = h(1).Faces(keep,:);

  
  
%% Show Array  
function arrange_axes_xy

figure(3), clf

plots.view_mesh('-f')

view([0 1 0]), axis([-3 3 -0.01 0.5 -3 3])

h = get(gca,'Children');

h(1).EdgeColor = [.3 .3 .3];
plot3(-[2 3],[0 0],[-3 -3],'k-','LineWidth',2,'clipping','off')

delete(h(3:end))

x = h(2).Vertices(:,3);
y0 = [ min(min(x(h(2).Faces)))  max(max(x(h(2).Faces))) ];

% C = lines(7);
% fill3([-3 3 3 -3 -3], [.2 .2 .2 .2 .2], y0([1 1 2 2 1]), C(1,:), 'FaceAlpha',0.5,'EdgeColor',C(1,:))
% delete(h(2))

% x = h(2).Vertices(:,1);
% keep = ~all( x(h(2).Faces) > 3.1 |  x(h(2).Faces) < -3.1 , 2); 
% h(2).Faces = h(2).Faces(keep,:);

x = h(1).Vertices(:,2);
keep = all( x(h(1).Faces) > -0.01 , 2); 
h(1).Faces = h(1).Faces(keep,:);

% h(1).Vertices( h(1).Vertices(:,2) < -0.01 ,2) = NaN;

% savefig_pdf array_mesh_PN.pdf
