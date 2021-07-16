


function preview_axons(filename,varargin)
% preview the axon arrangement in a group of fascicles 
% (extracted from models.axon_population)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if nargin == 0, filename = tools.file('get','axons~/ax*.mat'); end
if contains(filename,'~'), filename = tools.file(filename); end
if ~exist(filename,'file'), 
  filename = tools.file('axons~/axons*.mat','-prompt');  
end

if isempty(whos('pop','file',filename))
    u = load(filename,'nerve','axon_populations');
    pop = u.axon_populations;
    nerve = u.nerve;
    clear u
else load(filename,'nerve','pop');
end

if ~exist('nerve','var')
  nF = max(pop.unmyelinated_fascicle);
  nerve.coeffs = zeros(2,0,nF);
  nerve.fascicles = zeros(0,2,nF);
  nerve.outline = zeros(0,2,nF);
end

if ~isfield(nerve,'fascicles')  
  nerve.fascicles = nerve.outline;
  assert(~iscell(nerve.outline))
end

nF = size(nerve.coeffs,3);
ds_factor = 1; % no downsampling

if any(named('-down')), ds_factor = get_('-down'); end


%%
if any(named('-fig')), figure(get_('-fig')); 
else figure(1), p = get(gcf,'Position'); 
  if all(p(3:4) == [560 420]), set(gcf,'Position',p .* [0.8 1 1.3 1]), end
end

clf, hold on, % C = lines(nF);

labels = {'Fascicle outline','unmyelinated afferent','unmylinated efferent', ...
            'myelinated afferent','myelinated efferent'};

if max(nerve.fascicles(:)) > 20 % mm to um
  nerve.fascicles = nerve.fascicles / 1000; 
end         
if iscell(nerve.outline)
  fill(nerve.outline{end}(:,1),nerve.outline{end}(:,2),[.9 .9 .9],'EdgeColor','none')
  labels = [{''} labels];
end
for ff = 1:nF
    fill(nerve.fascicles(:,1,ff),nerve.fascicles(:,2,ff),'w','EdgeColor',[.3 .3 .3],'LineWidth',1.2)
    
    if ~any(named('-no-l')), 
      text(mean(nerve.fascicles(:,1,ff)),mean(nerve.fascicles(:,2,ff)),sprintf('\\bfF%d',ff), ...
          'Color',[.7 .7 .7],'FontSize',26,'Horiz','center')
    end

    for pp = numel(pop):-1:1
    
        f_id = (pop(pp).fascicle == ff); 
        xy = pop(pp).axon_xy;
    
        if pop(pp).myelinated, 
            style = {'o','MarkerSize',5, ...
                         'MarkerFaceColor',(pop(pp).color+1.5)/2.5', ...
                         'LineWidth',1.2};
        else style = {'.','MarkerSize',10};
        end
        
        f_id(f_id) = mod(1:sum(f_id), ds_factor) == 0;
        plot(xy(f_id,1),xy(f_id,2),style{:},'color',pop(pp).color)
    end
end

if any(named('-e')), e = get_('-e'); 
  if isfield(e,'info'), e = e.info; end
  
  y2x = range(ylim) / range(e.ElectrodePositions(:,1));
  
  for cc = 1:size(e.ElectrodePositions,1)
    
    xy = [e.ElectrodePositions(cc,[3 1]) 0 0] + [-1 -1 2 2]/2 .* ...
          e.ElectrodeDimensions(e.ElectrodeTypeIndex(cc),[3 1 3 1]);
    
    xy(2) = xy(2) - max(e.ElectrodePositions(:,1));
    xy([2 4]) = xy([2 4]) * y2x;
    
    rectangle('Position',xy,'FaceColor',[.6 .6 .6],'EdgeColor',[.4 .4 .4])
    
  end
end

axis image, tools.tidyPlot
title(regexprep(tools.file('T',filename),'([\\_\^])','\\$1'))
legend(labels{:},'location','best')

if any(named('-no')), return, end

%%

if any(named('-fig2')), figure(get_('-fig2')); 
else figure(2), p = get(gcf,'Position'); 
  if all(p(3:4) == [560 420]), set(gcf,'Position',[220 120 730 700]), end
end

%% Produce figure panels 
clf
xmax = round(quantile(cat(1,pop.fibre_diam),0.998)); 
bar_style = {'EdgeColor','none','FaceAlpha',0.5,'FaceColor'};

[~,x] = hist(cat(1,pop.fibre_diam),36); %#ok<HIST>

G = @(v) [v v v]/10;

for pp = 1:numel(pop)
  if pop(pp).myelinated
        subplot(3,3,[4 8]), hold on,
        plot(pop(pp).fibre_diam,pop(pp).axon_diam,'o','Color', ...
                                pop(pp).color,'MarkerSize',4)     
  end
end

%%

is_myel = find([pop.myelinated]);   
    
if any(is_myel)
    subplot(3,3,[4 8]),
    axis equal, tools.tidyPlot, axis([-0.2 xmax -0.2 xmax])
    h = flipud(get(gca,'children'));
    
    g = 0.1:0.1:1;
    plot([0;xmax;nan]*(0*g+1),[0;xmax;nan]*g,'-','Color',[0 0 0 0.3])
    text([xmax 0]*[1.01; -0.01],0.1*xmax,'g=0.1','Color',G(7))
    text([xmax 0]*[1.01; -0.01],0.9*xmax,'g=0.9','Color',G(7))
    grid on

    legend(h,pop(is_myel).axon_type,'location','nw','box','off')

    xlabel('fibre diameter (µm)'),
    ylabel('axon diameter (µm)'),
    
    
    subplot(3,3,[1 2]), hold on        
    y = hist(cat(1,pop(is_myel).fibre_diam),x); %#ok<HIST>
    bar(x,y,1,bar_style{:},pop(is_myel(end)).color)
    
    xlim([0 xmax]), tools.tidyPlot, ylabel('myelinated')
    set(gca,'YColor',pop(is_myel(end)).color)
    ylim([0 max(ylim)])

    if any(~[pop.myelinated]), axes('position',get(gca,'Position')), end
else cla
end

is_unmy = setdiff(1:pp,is_myel); 

y = hist(cat(1,pop(is_unmy).fibre_diam),x); %#ok<HIST>
bar(x,y,1,bar_style{:},pop(is_unmy(1)).color)
    
xlim([0 xmax]), tools.tidyPlot,  ylabel('unmyelinated')
set(gca,'YColor',pop(is_unmy(1)).color,'YAxisLocation','right')
ylim([0 max(ylim)])


if any(is_myel), subplot(3,3,[6 9]),

    [y,x] = hist(cat(1,pop(is_myel).g_ratio),0:0.05:1); %#ok<HIST>
    barh(x,y,1,bar_style{:},pop(is_myel(end)).color)
    tools.tidyPlot, ylabel('g-ratio')
    set(gca,'YTick',0:0.1:1)
    xlabel('axon count')
end

%%
return
