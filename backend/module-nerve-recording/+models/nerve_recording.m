
function settings = nerve_recording(varargin)
% TODO add header documentation : raw docstring below
%
%  if any(named('-list-modes')), return settings
%  if any(named('-file')), eidors_file = get_('-file'); end
% 
%  if ('-recenter-peak')) % explan this complicated proc
%  if ('-root')), opts.axons_folder = get_('-root'); end
%  if ('-dis')),
%  if ('-debug-distortion')), opts.plot_image = true; end
%  if ('-coh')), cohere_settings = get_('-coh');
%  if ('-zref')) z_ref = get_('-zref'); end
%  if ('-fs')),      fs = get_('-fs');
%  if ('-time')),  time_span = get_('-time');
%  if ('-rep')), settings.n_reps = get_('-rep'); end
%  if ~('-unit-um')), xy = xy / 1e3; end % um -> mm
%  if ('-set-r')), get_raster_ = get_('-set-r');
%  if ('-debug-u'))
%  if ('-DEBUG'))
%  if ('-fascicle-sum'))
%  if ('-out')), e_name = get_('-out');
% 
% -f-list [list] : simulate just the listed fascicles 
% 
% in function opts = get_settings(varargin):
%  named = @(v) strncmpi(v,varargin,length(v));
%  get_ = @(v) varargin{find(named(v))+1};
%  if any(named('-settings')), opts = get_('-settings'); return, end
%  if any(named(['-' opts(x).name])), opts = opts(x); return, end
% 
% opts(1).name = 'default';
% opts(2).name = 'base';
% opts(3).name = 'drift';
% opts(4).name = 'ecap';
% opts(5).name = 'burst';
% opts(6).name = 'phase';
% opts(7).name = 'flat';
% 

varargin = tools.opts_to_args(varargin,'recording');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('-q')), printf = @(varargin) []; else printf = @fprintf; end

if ~isdeployed && isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP>
    save_default_options ('-mat-binary'), 
end

settings = get_settings(varargin{:}); 

if any(named('-list')), return    
else settings = settings(1); 
end

printf('Running models.%s ... \n', mfilename);

[EM,AX] = tools.parse_arguments(varargin, 'LOAD', ...
                          'eidors','eidors~/sens*.mat', 'axons');
opts = struct;
opts.axons_folder = AX.folder; 
opts.verbose = ~any(named('-q'));
opts.use_parallel = ~any(named('-no-p'));

nerve = AX.nerve;
pop = AX.pop;

for f = fieldnames(EM.utils)' % Create utility local functions
  EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');
  eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
end

%% Get sensitivity functions 

fascicles = fieldnames(EM); % "Fascicle1" not guaranteed to exist
fascicles(~strncmp(fascicles,'Fascicle',8)) = [];

nE = size(EM.(fascicles{1}).pot,1);
nF = size(nerve.fascicles,3); 

sensitivity = cell(nE,nF);
for ff = 1:nF % for each fascicle and electrode, make i2v function 
  fascicle = sprintf('Fascicle%d',ff);     
  if ~isfield(EM,fascicle), continue, end
  for ee = 1:nE 
   if ee == 1, sensitivity{ee,ff} = scatteredInterpolant( ... 
                       z_(fac_(ff)), y_(fac_(ff)), x_(fac_(ff)), ...
                       EM.(fascicle).pot(ee,:)', ...
                       'linear','none'); % units of mm ? 
   else sensitivity{ee,ff} = sensitivity{1,ff}; 
        sensitivity{ee,ff}.Values = EM.(fascicle).pot(ee,:)';
   end
 end
end

clear ee ff x_ y_ z_ fac_ f sel ok

% Translate sensitivity peaks if specified
if any(named('-recenter-peak'))
  disp('Translating sensitivity peaks...')
  sensitivity = translate_I2V_peaks(sensitivity, get_('-recenter-peak'));
end

%% Get Distortion if specified

if any(named('-dis'))
  opts.distortion = get_('-dis'); 
  if size(opts.distortion,1) == numel(pop)
    opts.distorion_by_type = opts.distortion; 
  end
  % options.distortion. [ v z t ] is read as a matrix or structure by
  % models.spike_to_wave.
  
  if any(named('-debug-dis')), opts.plot_image = true; end  
end

%%

if strcmpi(settings.name,'ECAP'), z_ref = 'E3'; else z_ref = 0; end

if any(named('-coh')), cohere_settings = get_('-coh');
else  cohere_settings = settings.coherence;
end

population_frequency = settings.frequency;
population_spikerate = settings.spikerate;
population_exponent  = settings.exponent;

% default position for spiketime: a nominal spike-time of "0" 
%   means the AP is at z=0 at time 0 (axons oriented along
%   z-axis and propegate in the +z direction) 
if any(named('-ref')) z_ref = get_('-ref'); end
if ischar(z_ref), z_ref = str2double(regexp(z_ref,'\d+','match','once'));
    z_ref = median(EM.model.nodes(EM.model.electrode(z_ref).nodes,1)); 
end

if any(named('-axon-t')), EM.info.AxonTrajectory = get_('-axon-t'); end

fs = 30; % 24.414; % kHz 
if any(named('-sample-rate')), fs = get_('-sample-rate'); 
elseif any(named('-fs')),      fs = get_('-fs'); 
end

time_span = 65;
if any(named('-time')),  time_span = get_('-time'); 
 if strcmp(time_span,'auto'), time_span = max(results(:,3) ./ results(:,2)); end %#ok<NODEF>
end

if any(named('-rep')), settings.n_reps = get_('-rep'); end
n_rep = settings.n_reps;

if any(named('-raster')),
  input_raster = get_('-raster');
  if ~isa(input_raster,'function_handle')
  
    if ischar(input_raster), 
      printf('Loading %s ... \n', input_raster);
      input_raster = load_spikes_file(input_raster, fs); 
    end

    population_frequency = 1; 
    population_spikerate = 1:numel(input_raster); 
    population_exponent = 1; 
    n_rep = 1; 
    
    settings.name = 'input'; 
    settings.wave_path = 'file';
    settings.file_scheme = 'epoch_%03d (%%d).mat';
    settings.file_vector = @(ex,sr,fr,ch) sr(1);
    
  end
end
    
    

warn_once = true; 
check_folder = true; 

%% Compute intracellular-to-extracellular relationship
for i_rep = 1:n_rep
 for i_coh = 1:length(cohere_settings)
  for i_freq = 1:length(population_frequency)
   for i_rate = 1:size(population_spikerate,1)

    if size(population_spikerate,2) >= 2 && ...
            population_spikerate(i_rate,1) == ...
            population_spikerate(i_rate,2) && i_freq > 1
      continue
    end

    waves = []; 
    raster = {}; 
    options = []; 
    
    % load raster input here ????? 
    
    for ty = 1:numel(pop)
        
      time = 0:(1/fs):(time_span);
      time = [-fliplr(time) time(2:end)]; % ms

      opts.filename = '';      
      opts.reference_z = z_ref; 
      opts.class = pop(ty).axon_model;
      opts.padding = 2;
      opts.raster_opts.ph = 0; % 2*pi*rand      
      
      opts.raster_opts = struct; % defaults
      if strcmpi(settings.name,'drift')
        opts.raster_opts.ty = 'step';
        opts.raster_opts.ex = population_exponent(i_freq);
        opts.raster_opts.fc = 1;
      elseif strcmpi(settings.name,'burst')
        opts.raster_opts.ty = 'cos';
        opts.raster_opts.ex = population_exponent(1);
        opts.raster_opts.fc = population_frequency(i_freq);
      elseif strcmpi(settings.name,'phase')
        opts.raster_opts.ty = 'cos';
        opts.raster_opts.ex = population_exponent(1);
        opts.raster_opts.fc = population_frequency(i_freq);
        opts.raster_opts.ph = deg2rad(population_spikerate(i_rate,3));
      elseif strcmpi(settings.name,'pulse')
        opts.raster_opts.ty = 'exp2'; % double exponential
        opts.raster_opts.tau1 = 0.5 ;
        opts.raster_opts.tau2 = 2.5 ; % ms time constants
        opts.evoked_potential = true; % prevent currents at time < 0
      elseif strcmpi(settings.name,'gauss')
        if i_freq > 1, time = []; break, end
        opts.raster_opts.ty = 'gauss'; % gaussian spike-time profile
        opts.raster_opts.ex = population_exponent(1);
        opts.raster_opts.fc = population_frequency(i_rate);
      end
      
      if strcmpi(settings.name,'input')        
        opts.raster_opts = input_raster.epoch_info;
        if isfield(opts.raster_opts, 'z_ref')
          opts.reference_z = opts.raster_opts.z_ref; 
        end
      elseif size(population_spikerate,2) >= 2 
        opts.raster_opts.fb = population_spikerate(i_rate,1); % imp/s base 
        opts.raster_opts.fp = population_spikerate(i_rate,2); % imp/s peak  
      elseif strncmpi(settings.wave_path,'flat',4)
        opts.raster_opts.fb = population_spikerate(i_rate);
        opts.raster_opts.fp = population_spikerate(i_rate);
      else
        opts.raster_opts.fb = 0.1; 
        opts.raster_opts.fp = population_spikerate(i_rate);
      end
      
      opts.raster_opts.ax_sd = cohere_settings(i_coh); 
      
      xy = pop(ty).axon_xy;
      grp = pop(ty).size_sample;
      fid = pop(ty).fascicle;
      
      if any(named('-fasc')) % -fascicle-list
        ok = ismember(fid,get_('-fasc'));
        xy = xy(ok,:); grp = grp(ok); fid = fid(ok); 
      end
      
      if isfield(opts,'distortion_by_type')
        opts.distortion = opts.distortion_by_type(ty,:); 
      end
      
      % Work out whether the xy units need to be corrected
      if max(abs(xy(:))) > max(EM.model.nodes(:,1)),  xy = xy / 1e3;
          warning('ViNERS:axonUnits_mm',['%s should be in units of mm. ' ...
            'treating this as um, maybe cancel this and fix?'], ...
            tools.file('T',[opts.axons_folder opts.class '\index.mat']))
      elseif any(named('-unit-um')) xy = xy / 1e3;
      end % um -> mm ? 
    
      nG = max(grp);
      
      if isfield(settings,'modulation_fun')
        opts.raster_opts.modulation = settings.modulation_fun;
      end

      if strcmpi(settings.name,'input')
        if iscell(input_raster), opts.raster = input_raster{i_rate};
        else opts.raster = input_raster(i_rate);
        end
        
        pidx = [0 cumsum(cellfun(@numel,{pop(1:ty).fascicle}))]; 
        if isfield(opts.raster,'population_id') % identifying field
          sel = opts.raster.population_id == ty;
        elseif any(mod(opts.raster.spk_axon,1)) % id.pop# syntax
          sel = mod(opts.raster.spk_axon,1)*10 == ty;
        else          
          sel = opts.raster.spk_axon <= pidx(end) & ...
                opts.raster.spk_axon > pidx(end-1);
        end
        
        opts.raster.spk_axon = opts.raster.spk_axon(sel) - pidx(end-1);
        opts.raster.spk_time = opts.raster.spk_time(sel);
        
        if isempty(opts.raster.bin_rate), 
          opts.raster.bin_rate = opts.raster.pop_rate; 
        end
        
      elseif any(named('-raster')), 
        get_raster_ = get_('-raster');
        opts.raster_opts.loop_indices = [i_coh i_freq i_rate i_rep ty]; 
        opts.raster = get_raster_(time,xy,opts.raster_opts); 
      else opts.raster = models.random_raster(time,xy,opts.raster_opts);          
      end
      
      gff = 10.^ceil(log10(nF+2)); % subgroup_id = gg + ff / 10 would fail for nF > 9
      opts.raster.axon_group = grp + fid/gff; 
      opts.efferent = ~pop(ty).afferent; % propegation is in +Z direction for afferents, -Z for efferents 
      opts.spike_rate = mean(opts.raster.bin_rate);
      
      % plot(opts.raster.bin_time, opts.raster.bin_rate)
      % hold on, result = [result opts.spike_rate];

      if isfield(opts,'evoked_potential') && opts.evoked_potential 
        opts.efferent = false;
      end

      % axon_xy is a {nG x nF} cell array of axon positons. By default, 
      % each element is a [nCells x 2] list of intercepts in the xy plane.
      % If EM.info.AxonTrajectory or EM.info.FascicleTrajectory is defined,
      % each element is instead a [nCells x 3 x nPoints] trajectory for
      % each cell. 
      axon_xy = cell(nG,nF); 
      
      if tools.from_trajectory(EM,nerve)
        for ff = 1:nF % For each fascicle
          axon_xy(:,ff) = tools.from_trajectory(EM,nerve, ...
                                             xy(fid == ff,:),'-f',ff, ...
                                       '-g',grp(fid == ff),'-nG',nG);          
        end
      else
        for ii = 1:(nG*nF) 
          [gg,ff] = ind2sub([nG nF],ii); 
          axon_xy{ii} = xy(grp == gg & fid == ff,:); 
        end % setup g_xy    
      end

      V = cell(nG,nF); % Observed potentia; 
      % V1 = cell(nG,nF); ~ second argument is potential of first spike

      if any(named('-debug-u')) || any(named('-preview'))
        %% Debug this, is a consistent problem ... 

        clf, plots.view_mesh(EM.model,'-fasc')
        h = get(gca,'Children');
        h(2).EdgeAlpha = 0.25;  h(2).FaceAlpha = 0; 
        h(2).Clipping = 'off';
        
        for gg = 1:nG*nF
          for uu = 1:size(axon_xy{gg},1)
            xyz = permute(axon_xy{gg}(uu,:,:), [3 2 1]);
            if size(axon_xy{gg},2) == 2
              plot3(0*xyz(1,:), xyz(2,:),xyz(1,:),'o', ...
                  'Color',[1 0 0 0.2],'Clipping','off','markersize',5) 
            else
              plot3(xyz(:,1), xyz(:,2),xyz(:,3), ...
                  'Color',[1 0 0 0.2],'Clipping','off') 
            end
            hold on
          end
        end
        
        
        lim = [min(cat(1,axon_xy{:})) max(cat(1,axon_xy{:}))];
        if size(lim,3) > 1, lim = [min(lim(:,[3 2 1],:),[],3) ...
                                   max(lim(:,[6 5 4],:),[],3)];
        end
        if size(axon_xy{gg},2) > 2
          lim = [min(lim(:,1:2,:),[],3) max(lim(:,4:5,:),[],3)];
        end
        view([1 0 0]), 
        ylim(lim([2 4])*[1.1 -.1; -.1 1.1])
        zlim(lim([1 3])*[1.1 -.1; -.1 1.1])
        
        %%
        return
      end
      
      %% CORE parfor: Compute summation of currents from spike raster
      
      printf('Computing spatial summation of %d APs ... ',numel(opts.raster.spk_time))
      cache_path = tools.cache('path');

      if nG*nF == 1 || ~opts.use_parallel
        for ii = 1:(nG*nF)
          V{ii} = parfun_unpack(cache_path, @(g,f) ...
                        models.spike_to_wave(g+f/gff,time,        ...
                                             sensitivity(:,ff),  ...
                                             axon_xy{g,f},opts), ...
                                          [nG nF],ii);
        end
      elseif isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP> % octave parallel
        if isempty(which('pararrayfun')), pkg load parallel; end 
        V = pararrayfun(nproc-1, @(a) parfun_unpack(cache_path, @(g,f) ...
                                         models.spike_to_wave(g+f/gff,time, ...
                                                        sensitivity(:,ff), ...
                                                       axon_xy{g,f},opts), ...
                                    [nG nF],a), 1:(nG*nF), 'Unif', false);
        error TODO_validate_this_for_octave
      else    
        gff = gff;  %#ok<ASGSL>
        parfor ii = 1:(nG*nF)  
          V{ii} = parfun_unpack(cache_path, @(g,f) ...
                        models.spike_to_wave(g+f/gff,time,        ...
                                             sensitivity(:,ff),  ...
                                             axon_xy{g,f},opts), ...
                                          [nG nF],ii); %#ok<PFBNS>
        end
      end
      printf('Done!\n')
        
      for ff = 1:nF        
        waves(:,:,ff,ty) = sum(cat(3,V{:,ff}),3); %#ok<AGROW>
        % TODO convert axon_wave based on partial model
      end
      raster{ty} = opts.raster;  %#ok<AGROW>
      
      if isempty(options), options = opts; else options(ty) = opts; end %#ok<AGROW>
    end % ty [1..4]

    if isempty(time), continue, end

    % This check prevents from having to open another instance
    summary = [sqrt(nanmean(waves(:).^2)) nanmin(waves(:)) nanmax(waves(:))];
    printf(' computed wave: RMS %0.1f, range %01f-%0.1f\n', summary ); 
    if warn_once && all(summary == 0)
        warning('ViNERS:possibleUnitsFail','WAVE was all-zeros, try calling with -unit-um or -debug-units')
        warn_once = false; 
    end      
     
    
    %% Debug visualisations / inspection utilities 
    if 0
      %% Look at resulting overall recordings 

      clf %#ok<UNRCH>
      plot(time,waves(:,:,end) + (0:3)), ax = gca;
      xlim(time([1 end]))
      ax.Position(4) = ax.Position(4)/2;
      ax(2) = axes('Position',ax.Position +[0 1 0 0]*ax.Position(4));

      plot(raster{end}.spk_time,raster{end}.spk_axon,'k.')
      ax(2).XLim = ax(1).XLim;

      %% Look at outputs for individual spikes

      V_e1 = cellfun(@(v)v(:,4),vSpk,'UniformOutput',false);
      V_e1 = [V_e1{:}];

      spk1_t = nan*V_e1(1,:);

      for ii = 1:(nG*nF) % get first spiketime 

       [gg,ff] = ind2sub([nG nF],ii); 
        g_id = gg+ff/10;
        ax_g = (opts.raster.axon_group == g_id); % Group mask
        ax_k = ax_g(opts.raster.spk_axon); % which times are in group?
        if any(ax_k), ax_k = find(ax_k,1);      
          spk1_t(ii) = opts.raster.spk_time(ax_k); 
        end
      end

      clf, hold on
      plot(time,V_e1 + (1:size(V_e1,2))/2,'Color',[0 0 0 0.5])
      plot(spk1_t,(1:ii)/2,'r.')
      tools.tidyPlot, axis tight, ylim([-0.5 max(ylim)])
      set(gca,'YTick',(0:nG:(nG*nF))/2)

      clear V_e1 gg ff g_id ax_g ax_k spk1_t
      % [~,ii] = ginput(1); ii = round(ii*2);

      %% I've gotten the impression something funny is happening? 

      opts.plot_wave = 1; 
      [gg,ff] = ind2sub([nG nF],ii); 
      [V_check,vSpk_check] = models.spike_to_wave(gg+ff/10, ... 
                                    time,sensitivity(:,ff),axon_xy{ii},opts);

      %% Check the saved "results" field
      nG = numel(sam.A_diam);
      for gg = 1:nG % Load cache results 
        d = load(axonpath(ty,'n%03d_out.mat',gg),'result');      
        results(gg,:) = d.result;
      end

    end % if plot_debug

   %%
    if any(named('-f-sum')), waves = squeeze(sum(waves,3)); 
    elseif size(waves,3) == 1, waves = squeeze(waves); 
    end

    clear gg ff ee g_xy ii xy grp ty fid ans
    
    if any(named('-out')), wave_path = get_('-out');
      if any(wave_path == '~'), wave_path = tools.file(wave_path); end
      if ~any(ismember('\/',wave_path)) % put in ~waves/
        if any(ismember('()',wave_path)) % cat settings.wave_path 
          wave_path = [settings.wave_path ' ' wave_path];  %#ok<AGROW>
        end
        wave_path = [tools.file('waves~/') wave_path ]; %#ok<AGROW>
      end
    else % replace sensitivity with wave_path
      wave_path = regexprep(EM.filename,'sens[^\s\(]*', settings.wave_path);
      wave_path = [tools.file('waves~/') wave_path ]; %#ok<AGROW>
    end
    
    wave_path = [wave_path filesep]; %#ok<AGROW>
    
    file_out = sprintf(settings.file_scheme,settings.file_vector( ...
                      settings.exponent(min(i_freq,end)), ...
                      settings.spikerate(i_rate,:), ...
                      settings.frequency(min(i_freq,end)), ...
                      settings.coherence(i_coh))); 
    
    file_out = tools.file('get',[wave_path strrep(file_out,'.0_','_')],'next'); 

    if ~exist(fileparts(file_out),'dir'), mkdir(fileparts(file_out)),
      check_folder = false; 
    elseif check_folder
      
      if any(named('-clf'))
          warning('ViNERS:overwriteFolder','erasing the pre-existing directory %s', tools.file('T',wave_path))
          rmdir(fileparts(file_out),'s'), mkdir(fileparts(file_out))
      else
          warning('ViNERS:possibleOverWrite','possibly writing in a pre-existing directory %s', tools.file('T',wave_path))
      end
      check_folder = false; 
    end
    
    
    
    options(1).class = {options.class};
    options(1).spike_rate = [options.spike_rate];
    options = rmfield(options(1),{'verbose','raster','efferent'});
    options.afferent = [pop.afferent];
    options.filename = file_out;  
    inputs = varargin;    
    raster = [raster{:}]; 
    
    if isdeployed,
      n_total = n_rep * length(cohere_settings) *  size(population_spikerate,1);
      n_done = numel(dir([wave_path '/*.mat']));         
      fprintf('Progress: %0.1f%%\n', 10+75*(n_done+1)/n_total);
    end

    
    % printf('Saving %s\n', tools.file('TT',file_out))
    printf('Saving %s\n', tools.file('T',file_out))
    save(file_out, 'raster','waves','time','inputs','options')

   end % i_rate
  end % i_freq
 end % do_population
end % do_replicates


if nargout == 0, clear, end
return


function opts = get_settings(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

opts = struct; 

if any(named('-settings')), opts = get_('-settings'); return, end

opts(1).name = 'flat';
opts(1).n_reps = 3;
opts(1).coherence = 1;
opts(1).spikerate = [0.1; 1; 10];
opts(1).frequency = 1;
opts(1).exponent = 3; 
opts(1).wave_path = 'flat';
opts(1).file_scheme = 'epoch_k%0.1f_c%0.1f (%%d).mat';
opts(1).file_vector = @(ex,sr,fr,ch) [sr(1) ch]; 

opts(2) = opts(1); 
opts(2).name = 'drift'; 
opts(2).frequency = [1 1 1];
opts(2).spikerate = [0.1 1; 1 5; 1 10; 1 20; 3 4];
opts(2).exponent  = [1 5 20];  
opts(2).coherence = [0.3 1 3]; 
opts(2).wave_path = 'drift'; 
opts(2).file_scheme = 'epoch_b%0.1f_k%0.1f_c%0.1f_w%0.0f (%%d).mat';
opts(2).file_vector = @(ex,sr,fr,ch) [sr(1:2) ch ex]; 

opts(3) = opts(1); 
opts(3).name = 'burst';
opts(3).frequency = [100 50 20 10];
opts(3).spikerate =  [.1 3; 0.75 1.5; 1 1; 1 30; 7.5 15; 10 10 ]; 
opts(3).wave_path = 'burst'; 
opts(3).file_scheme = 'epoch_k%0.1f_f%0.1f_c%0.1f (%%d).mat';
opts(3).file_vector = @(ex,sr,fr,ch) [sr(2) fr ch]; 

for x = 1:numel(opts)  
  if any(named(opts(x).name)), opts = opts(x); return, end  
  if any(named(['-' opts(x).name])), opts = opts(x); return, end  
end

return

%% Move sensitivity peak (for electrode array property explorations) 
function sensitivity = translate_I2V_peaks(sensitivity, z_REF)
%%

if numel(z_REF) < nE
  z_REF = repmat(z_REF(:), [ceil(nE/numel(z_REF)) 1]); 
end

color = 1-summer(nE);
figure(3), clf

for ff = 1:nF % for each fascicle 
  xy = [mean(z_(fac_(ff))) mean(y_(fac_(ff)))];
  z = linspace(min(x_(fac_(ff))), max(x_(fac_(ff))),1001);
  subplot(nF,1,ff), hold on
  for ee = 1:nE

    pk0 = sensitivity{ee,ff}([xy 0] + z'*[0 0 1]); 
    [~,id] = nanmax(pk0);   
    sensitivity{ee,ff}.Points(:,3) = sensitivity{ee,ff}.Points(:,3) ...
                                     - z(id) + z_REF(ee);
    pk1 = sensitivity{ee,ff}([xy 0] + z'*[0 0 1]);    
    plot(z,pk0,':','Color',color(ee,:),'LineWidth',1.2)
    plot(z,pk1,'-','Color',color(ee,:),'LineWidth',1.2)
    % pause(0.02)
  end
  tools.tidyPlot, pause(0.05)
end 

return

function wave = parfun_unpack(cache_path,fun,nGnF,ii)
  tools.cache('set',cache_path)
  [gg,ff] = ind2sub(nGnF,ii);
  wave = fun(gg,ff); 

function raster = load_spikes_file(filename, fs)

has_ext_ = @(b) strncmpi(fliplr(filename),fliplr(b),length(b)); 
if has_ext_('.mat'), raster = load(filename);
elseif has_ext_('.xml'), raster = convert_xml_raster(filename, fs); 
elseif has_ext_('.json'), raster = tools.parse_json(filename);
elseif has_ext_('.nwb'), raster = tools.read_NWB_file(filename);
else error('Unknown file format "%s", expected {.mat,.json,.xml}');
end
  
raster = validate_raster(raster, fs);

return

function raster = validate_raster(d, fs)

if numel(d) > 1, 
  if iscell(d), raster = cellfun(@validate_raster,d);
  else raster = arrayfun(@validate_raster,d);
  end, return
end

spk_unit = []; 
for u = {'spk_axon','axon','unit','u','cell','ax','id'}
  if isfield(d,u{1}), spk_unit = d.(u{1}); break, end
end

spk_time = []; 
for u = {'spk_time','time','t'}
  if isfield(d,u{1}), spk_time = d.(u{1}); break, end
end

info = struct;
for u = {'epoch_info','file_info','info','metadata','meta'}
  if isfield(d,u{1}), info = d.(u{1}); break, end
end

assert(~isempty(spk_time),'missing spike-times');

nA = []; 
if isfield(info,'n_units'), nA = info.n_units;
elseif isfield(info,'n_axons'), nA = info.n_axons;
elseif isfield(info,'n_cells'), nA = info.n_cells;
elseif isfield(info,'n_spikes'), nA = info.n_spikes;
elseif ~isempty(spk_unit), nA = nanmax(spk_unit); 
end
assert(~isempty(nA), 'missing n_axons'); 

if isempty(spk_unit), 
  disp('Assigning random spike units');
  spk_unit = ceil(rand(size(spk_time))*nA); 
end

if nargin < 2, fs = 30; end
for u = {'output_sample_rate','sample_rate','fs'}
  if isfield(info,u{1}), fs = info.(u{1}); break, end
end

ok = ~isnan(spk_time) & ~isnan(spk_unit); 
if ~all(ok)
  warning('%d/%d spiketimes corrupted.',sum(~ok), numel(ok))
  spk_time = spk_time(ok); 
  spk_unit = spk_unit(ok);
end


if isfield(info,'time_start') && isfield(info,'time_end')  
     time = (info.time_start) : (1/fs) : (info.time_end);   
else time = floor(min(spk_time)) : (1/fs) : ceil(max(spk_time)); 
end

if isfield(info,'bin_dt')    
     bin_x = min(time):(info.bin_dt):max(time);
else bin_x = min(time):(1):max(time); % default 1 ms bins
end, bin_y = hist(spk_time,bin_x);  %#ok<HIST>
bin_y = bin_y ./ mean(diff(bin_x)) / nA * 1000;

raster.epoch_info = info; % standardise variable names 
raster.spk_time = spk_time;
raster.spk_axon = spk_unit;

raster.spk_rate = bin_y;
raster.bin_time = bin_x; 
raster.bin_rate = []; % undefined
raster.pop_rate = numel(spk_time) / range(time) / nA * 1000;

return

function raster = convert_xml_raster(filename, fs)

if ~exist('filename','var'), filename = 'example.xml'; end
xml = tools.parse_xml(filename);  
xml.filename = filename;

%% x.Children().name representation
the = @(x,n) x.Children(strncmpi({x.Children.Name},n,length(n)));
attr_ = @(x,n) x.Attributes(strcmpi(n,{x.Attributes.Name})).Value;
trim_ = @(x) x.Children(~contains({x.Children.Name},'#text'));

xml.Children = trim_(xml); 

info = the(xml,'config'); 
info.Children = trim_(info);

info.Data = struct; 
info.Data.epoch_id = 0; 
info.Data.filename = filename;
info.Data.syntax = xml.Name;

for ii = 1:numel(info.Children)
  val = attr_(info.Children(ii),'value'); 
  if ~isnan(str2double(val)), val = str2double(val); end
  info.Data.(attr_(info.Children(ii),'name')) = val; 
end, shared_info = info.Data; 

shared_info.output_sample_rate = fs;

epochs = the(xml,'epoch');

if isempty(epochs), epochs = xml; end

[epochs.Data] = deal([]); 

%% Get spike-times

for ee = 1:numel(epochs)

  info = shared_info;
  info.epoch_id = ee; 

  for ii = 1:numel(epochs(ee).Attributes)
    val = epochs(ee).Attributes(ii).Value; 
    if ~isnan(str2double(val)), val = str2double(val); end
    info.(epochs(ee).Attributes(ii).Name) = val;
  end

  pop_block = the(epochs(ee),'pop'); 
  if isempty(pop_block)        
    epochs(ee).Data = parseXMLPopulation(epochs(ee),info,[]);
  else
    clear pop_data
    for pp = 1:numel(pop_block)
      pop_data(pp) = parseXMLPopulation(pop_block(pp),info,pp); %#ok<AGROW>
    end
    for u = fieldnames(pop_block(pp).Data)'
      epochs(ee).Data.(u{1}) = cat(1,pop_data.(u{1})); 
    end
    epochs(ee).Data.epoch_info = epochs(ee).Data.epoch_info(1);
  end
end

raster = cat(1,epochs.Data); 

return

function raster = parseXMLPopulation(xml,info,pop_id)

the = @(x,n) x.Children(strncmpi({x.Children.Name},n,length(n)));
has_ = @(x,n) any(strcmpi(n,{x.Attributes.Name}));
attr_ = @(x,n) x.Attributes(strcmpi(n,{x.Attributes.Name})).Value;

spk_time = []; 
spk_unit = []; 
spk_pidx = []; 

spikes = the(xml,'s');

for ii = 1:numel(spikes)
  spk_time = [spk_time; str2double(attr_(spikes(ii),'t'))]; %#ok<AGROW>
  spk_unit = [spk_unit; str2double(attr_(spikes(ii),'u'))]; %#ok<AGROW>
  if ~isempty(pop_id), spk_pidx = [spk_pidx; pop_id]; continue, end %#ok<AGROW>
  if ~has_(spikes(ii),'p'), continue, end
  spk_pidx = [spk_pidx; str2double(attr_(spikes(ii),'p'))]; %#ok<AGROW>      
end

units = the(xml,'unit');

for ii = 1:numel(units)
  uid = str2double(attr_(units(ii),'id'));
  if ~isempty(pop_id), pid = pop_id; 
  elseif has_(units(ii),'type'),
    pid = str2double(attr_(units(ii),'type'));
  elseif has_(units(ii),'axon_type'),
    pid = str2double(attr_(units(ii),'axon_type'));
  elseif has_(units(ii),'pop'),
    pid = str2double(attr_(units(ii),'pop'));
  else pid = []; 
  end
  t = textscan(units(ii).Children(1).Data,'%f');

  spk_time = [spk_time; t{1}];       %#ok<AGROW>
  spk_unit = [spk_unit; 0*t{1}+uid]; %#ok<AGROW>
  if isempty(pid), continue, end
  spk_pidx = [spk_pidx; 0*t{1}+pid]; %#ok<AGROW>
end

raster.epoch_info = info;
raster.spk_time = spk_time;
raster.spk_axon = spk_unit;
if ~isempty(spk_pidx), raster.population_id = spk_pidx; end

epochs.Data = raster; 
    % epochs.Data = validate_raster(raster); 