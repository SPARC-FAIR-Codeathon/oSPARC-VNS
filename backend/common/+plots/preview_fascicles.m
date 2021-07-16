

function s = preview_fascicles (varargin)
% preview display for mesh.insert_gmsh_fascicles
% TODO: add visualisation for 3d fascicles

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

s = mesh.insert_gmsh_fascicles(varargin{:},'-info','-all');

if isfield(s,'splines'), s = s.splines; end

if ~any(named('-cla')), cla reset, hold on, end

if numel(s) > 0
    
  for ii = numel(s):-1:2
    plot_layer(s(ii),named,(1/ii).^0.2)
  end
  do_warn = plot_layer(s(1),named,1);
          

elseif iscell(s.coeffs)
  error('cell coeffs')
  if iscell(s.coeffs{1}), s.coeffs = s.coeffs{1}{varargin{1}};
                          s.outline = s.outline{1}{varargin{1}};
  else s.coeffs = s.coeffs{varargin{1}};
       s.outline = s.outline{varargin{1}};
  end
else
  do_warn = plot_layer(s,named);
end

axis equal, tools.tidyPlot, grid on, ax = axis; 

G = @(v) [v v v]/10;

if any(named('-array')), do_warn = []; 
  try array = get_('-array'); catch array = struct; end
  if isfield(array,'carrier'), axis(axis)      
      
    IDx = 0.55; IDy = 0.2; IDr = 0; % Defaults from cuff.geo.variable
    if isfield(array.carrier,'cuff_IDy'), IDy = array.carrier.cuff_IDy; end
    if isfield(array.carrier,'cuff_IDx'), IDx = array.carrier.cuff_IDx; end
    if isfield(array.carrier,'cuff_IDr'), IDr = array.carrier.cuff_IDr; end
    
    IDx = IDx * 1e3; IDy = IDy*1e3; IDr = IDr*1e3;
    
    x = cos(linspace(0,pi/2,31))*IDr;
    y = sin(linspace(0,pi/2,31))*IDr;
    f_ = @(u) fliplr(u); 
    
    x = [x+IDx/2-IDr f_(IDr-IDx/2-x) IDr-IDx/2-x f_(x+IDx/2-IDr)];
    y = [y+IDy/2-IDr f_(y+IDy/2-IDr) IDr-y-IDy/2 f_(IDr-y-IDy/2)];
    plot(x,y,'k-','LineWidth',1.4) 
    
  end  
else
  fill(ax([1 2 2 1 1]),[0 0 ax([3 3]) 0],G(3),'FaceAlpha',0.3,'EdgeColor',G(3),'LineWidth',1.1)
end

if any(named('-no-w')), return, end

if ~isempty(do_warn)
    do_warn = sprintf(',#%d',do_warn);
    text(mean(xlim),mean(ylim),sprintf('Negative Y in fascicles: %s', ...
                               do_warn(2:end)),'Color','r','FontSize',14)
end

if nargout == 0, clear, end
%%

function do_warn = plot_layer(s,named,darkness)

nF = size(s.coeffs,3);
C = lines(max(nF,7)); G = @(v) [v v v]/10;
do_warn = []; 

for ii = 1:nF
    
  if nargin < 3 || darkness == 1, c = C(ii,:);
  else c = G(min(3,5*darkness)); 
  end
  
  if any(named('-coe'))
    plot(s.coeffs(1,:,ii), s.coeffs(2,:,ii),'o','Color',c,'MarkerSize',4)
  end
  
  xy = s.outline(:,:,ii);
  if nargin > 2       
       fill(xy(:,1), xy(:,2),G(darkness*10),'EdgeColor',c,'LineWidth',1.2)  
  else plot(xy(:,1), xy(:,2),'-','Color',c,'LineWidth',1.2)
  end
  if any(s.outline(:,2,ii) < 0), do_warn = [do_warn ii]; end %#ok<AGROW>
  
end



