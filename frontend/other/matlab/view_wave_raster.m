
function view_wave_raster(varargin)
% This function generates the invididual panels for EMBC fig. 3

% TODO - documentation needs to go here. 

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

%% Figure illustrating responses [Single spike-rate example]

if nargin > 0, data_file = varargin{1}; else data_file = ''; end
if contains(data_file,'~'), data_file = tools.file(data_file); end
if ~exist(data_file,'file')
    data_file = tools.file('sub~/waves/*.mat','-prompt');
    if ~exist(data_file,'file'), return, end % cancelled
end

[data_dir,data_file] = fileparts(data_file); 
data_dir = [data_dir filesep];

D = load([data_dir data_file '.mat']);

ax_file = tools.file('get','sub~/axons/axon*.mat','newest');
if any(named('-ax')), ax_file = get_('-ax'); end
load(ax_file,'pop');

D.axontype = D.options.class;
D.axon_color = cat(1,pop.color);
% TODO - set electrode grid, plot opts for raster, types

if iscell(D.raster), D.raster = [D.raster{:}]; end

%%

t_roi = max(abs(D.time)) - 25; 
if any(named('-roi')), t_roi = get_('-roi');
  if ischar(t_roi), t_roi = str2double(t_roi); end
  if all(isnan(t_roi)), t_roi = 25; end
end
if numel(t_roi) == 1, t_roi = [-1 1]*abs(t_roi); end

clf, subplot(3,1,1), cla, hold on

z0 = 0; 
for ii = 1:numel(D.axontype)
  
    if isstruct(D.raster) && isfield(D.raster,'spike') % from models.ecap_recording
      spk_time = {D.raster(ii).spike.init}; 
      spk_axon = arrayfun(@(n) n*ones(size(spk_time{n})), 1:numel(spk_time),'unif',0);      
      spk_time = cat(1,spk_time{:}); spk_axon = cat(1,spk_axon{:}); 
      if isempty(spk_time), continue, end,       
      spk_time(:,2) = []; spk_axon(:,2) = []; 
    else
      spk_time = D.raster{ii}.spk_time;
      spk_axon = D.raster{ii}.spk_axon;
    end
        
    ok = (spk_time >= t_roi(1) & spk_time < t_roi(end));
    plot(spk_time(ok), spk_axon(ok) + z0,'.','Color', D.axon_color(ii,:))
    z0 = z0 + length(D.raster(ii).axon_group);       
end

tools.tidyPlot, set(gca,'YTick',[],'TickLength',[1 1]/150), xlim(t_roi)

% Row 2 - response waves

dy = quantile(abs(D.waves(abs(D.waves)  > 1e-12)),0.99); % µv

if any(named('-dy')), dy = get_('-dy'); end

if size(D.waves,4) > 1
   
  if any(named('-f')), fid = get_('-f'); else fid = 1:size(D.waves,3); end
  D.waves = permute(sum(D.waves(:,:,fid,:),3),[1 2 4 3]);
    
end


subplot(3,1,2), cla, hold on
style = {'linewidth',1.1,'color'};
C = lines(7); G = @(v) [v v v]/10; %#ok<NASGU>

ok = (D.time > t_roi(1)& D.time < t_roi(end));
plot(D.time(ok),dy+sum(D.waves(ok,3,:) - D.waves(ok,4,:),3),style{:},G(5))
plot(D.time(ok),sum(D.waves(ok,1,:) - D.waves(ok,2,:),3),style{:},G(2))

for ii = 1:numel(D.axontype)
  plot(D.time(ok),D.waves(ok,1,ii)-D.waves(ok,2,ii)-dy*ii,style{:},D.axon_color(ii,:))
end

axis tight, tools.tidyPlot, xlim(t_roi)
set(gca,'YColor','none','TickLength',[1 1]/150)

sb = min(10.^round(log10(max(-ylim))), 10);

plot((xlim*[1.01;-0.01])*[1 1],[0 sb]+min(ylim),style{:},get(gca,'XColor'),'Clipping','off')
text(xlim*[1.025;-0.025],sb/2+min(ylim),sprintf('%g µV', sb),'Color',get(gca,'XColor'),'Rotation',90, ...
    'horizontalalignment','center','verticalalignment','bottom')

xlabel('time (ms)')  
linkaxes(get(gcf,'children'),'x')
  
%%
subplot(3,1,3), cla reset, hold on

chronux_opts = tools.setupChronux; 
chronux_opts.Fs = 1000 / mean(diff(D.time)); % units of kS/s
chronux_opts.fpass = [0.2 4e3];

avg_spec = []; 
data_list = dir([data_dir regexprep(data_file,'k\d+','*')]);

if isempty(data_list)  
  data_list = dir([data_dir regexprep(data_file,'\(\d+\)','*')]);
end

if numel(data_list) <= 1
  warning('ViNERS:viewWaves:cannotAverage',... 
        'Only %d files found which are replicates of "%s", cannot produce an average.', ... 
        numel(data_list), data_file)
  if isempty(which('mtspectrumc'))
      
      [avg_spec,hz] = pwelch([squeeze(D.waves(:,1,:)-D.waves(:,2,:)) ...  
                              sum(D.waves(:,1,:)-D.waves(:,2,:),3)], ...
                               [],[],[],chronux_opts.Fs);
      avg_spec(hz > hz(end)/2,:) = []; 
       hz(hz > hz(end)/2,:) = []; 
      
      
  else [avg_spec,hz] = mtspectrumc([squeeze(D.waves(:,1,:)-D.waves(:,2,:)) ...  
                                        sum(D.waves(:,1,:)-D.waves(:,2,:),3)], ...
                                        chronux_opts);
  end
end

%%
for ff = 1:length(data_list)  
  
  D = load([data_dir data_list(ff).name]);
 
  D.axontype = D.options.class;
  D.axon_color = cat(1,pop.color);
  

  if size(D.waves,4) > 1

    if any(named('-f')), fid = get_('-f'); else fid = 1:size(D.waves,3); end
    D.waves = permute(sum(D.waves(:,:,fid,:),3),[1 2 4 3]);

  end

  if strcmp(data_list(ff).name,data_file), alpha = 0.6;
  else                                     alpha = 0.1;
  end
  
  if isempty(which('mtspectrumc'))
      [spec,hz] = pwelch([squeeze(D.waves(:,1,:)-D.waves(:,2,:)) ...  
                              sum(D.waves(:,1,:)-D.waves(:,2,:),3)], ...
                               [],[],[],chronux_opts.Fs);
      spec(hz > hz(end)/2,:) = []; 
       hz(hz > hz(end)/2,:) = []; 
  else [spec,hz] = mtspectrumc([squeeze(D.waves(:,1,:)-D.waves(:,2,:)) ...  
                                    sum(D.waves(:,1,:)-D.waves(:,2,:),3)], ...
                               chronux_opts);
  end

  for ii = 1:numel(D.axontype)
    plot(hz,spec(:,ii),'color',[D.axon_color(ii,:) alpha])
  end
  % plot(hz,spec(:,end),style{:},G(2))

  if isempty(avg_spec), avg_spec = spec;
  else avg_spec = avg_spec + spec;
  end
end

if ~isempty(data_list), avg_spec = avg_spec / ff; end

for ii = 1:numel(D.axontype)
  plot(hz,avg_spec(:,ii),style{:},D.axon_color(ii,:))
end
plot(hz,avg_spec(:,end),style{:},G(2))

xlabel('frequency, Hz'), ylabel('Spectral power')
set(gca,'YScale','log'), tools.tidyPlot

% ylim(10.^[-9.5 0])

tools.suptitle(strrep(data_file,'_','\_'))

  
%%
