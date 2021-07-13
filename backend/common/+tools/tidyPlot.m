% tidyPlotForIllustrator
% Sets up the current axes with pleasant looking values for publication
% graphs. These are intended for use with four subplots per page.

% Update PRM 
% Update CDE - can now take a handle to an axis as argument

function tidyPlot(options)

if nargin > 0 && isa(options(1),'matlab.graphics.axis.Axes')
     handle = options; options = struct;
else handle = gca;
end

if nargin < 1, options = struct; end

if ~isfield(options,'XminorTick'), options.XMinorTick = 'off'; end
if ~isfield(options,'Ygrid'),      options.Ygrid = 'off';      end
if ~isfield(options,'FigColor'),   options.FigColor = 'w';      end
    
  set(gcf,'Color',options.FigColor)

  set(handle, 'Box', 'off', 'Color', 'none',  ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.015 .015] , ...
              'YGrid'       , options.Ygrid, ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'LineWidth'   , 1.1);
             

%               'XMinorTick'  , options.XMinorTick, ...
          
          
if numel(handle) > 1, return, end

xl = get(handle,'xlim'); yl = get(handle,'ylim');

if strcmp(get(handle,'Xscale'), 'linear') && (xl(1) - xl(2)/30 < xl(2))
     set(handle, 'xlim', [xl(1) - xl(2)/30, xl(2)]);
end
if strcmp(get(handle,'Yscale'), 'linear') && (yl(1) - yl(2)/20 < yl(2))
       set(handle,'YLim', [yl(1) - yl(2)/20, yl(2)]); 
end
end
