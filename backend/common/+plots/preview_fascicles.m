

function s = preview_fascicles (varargin)
% preview display for mesh.insert_gmsh_fascicles

s = mesh.insert_gmsh_fascicles(varargin{:},'-info');

if isfield(s,'splines'), s = s.splines; end

if iscell(s.coeffs)
  if iscell(s.coeffs{1}), s.coeffs = s.coeffs{1}{varargin{1}};
                          s.outline = s.outline{1}{varargin{1}};
  else s.coeffs = s.coeffs{varargin{1}};
       s.outline = s.outline{varargin{1}};
  end
end

%%
cla reset, hold on

nF = size(s.coeffs,3);
C = lines(max(nF,7)); G = @(v) [v v v]/10;

do_warn = []; 

for ii = 1:nF
    plot(s.coeffs(1,:,ii), s.coeffs(2,:,ii),'o','Color',C(ii,:),'MarkerSize',4)
    plot(s.outline(:,1,ii), s.outline(:,2,ii),'-','Color',C(ii,:))   
    if nF > 7, 
        text(mean(s.outline(2:end,1,ii)), ...
             mean(s.outline(2:end,2,ii)),sprintf('#%d',ii),'Color',C(ii,:))
    end
    
    if any(s.outline(:,2,ii) < 0), do_warn = [do_warn ii]; end %#ok<AGROW>
end

axis equal, grid on, tools.tidyPlot, ax = axis; 

fill(ax([1 2 2 1 1]),[0 0 ax([3 3]) 0],G(3),'FaceAlpha',0.3,'EdgeColor',G(3),'LineWidth',1.1)

if ~isempty(do_warn)
    
    do_warn = sprintf(',#%d',do_warn);
    text(mean(xlim),mean(ylim),sprintf('Negative Y in fascicles: %s', ...
               do_warn(2:end)),'Color','r','FontSize',14,'Horiz','center')
end

if nargout == 0, clear, end