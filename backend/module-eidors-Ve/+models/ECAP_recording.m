
function ECAP_recording(varargin)
% TODO add header documentation : raw docstring below
% models.ECAP_response( sensitivity(*).mat, '-stim' stim_folder, ... )
% 
%  if any(named('-file')), eidors_file = get_('-file'); end
%  if ('-root')), opts.axons_folder = get_('-root'); end
%  if ('-root')), opts.stim_folder = get_('-stim'); end
%  if ('-out')), e_name = get_('-out');
% 
%  if ('-distortion')), apply systemic distortion to Im [v z t]
%  if ('-debug-distortion')), opts.plot_image = true; end
%  if ('-recenter-peak')) % explan this complicated proc
% 
%  if ('-fs')),      fs = get_('-fs');
%  if ('-time')),  time_span = get_('-time');
%  if ~('-unit-um')), xy = xy / 1e3; end % um -> mm
%  if any(named('-preview'))
%  if ('-fascicle-sum'))

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};
p_ = @(x) [x.folder filesep x.name];

if any(named('-q')), printf = @(varargin) []; else printf = @fprintf; end

[EM,AX,stim_folder] = tools.parse_arguments(varargin, 'LOAD', ...
                          'eidors','eidors~/sens*.mat', 'axons', 'stimulation');
options = struct;
% if any(named('-use-opts')), opts = get_('-use-opts'); end
options.axons_folder = AX.folder; 
options.stim_folder  = stim_folder;
options.verbose = ~any(named('-q'));
options.use_parallel = ~any(named('-no-p'));

printf('Running models.%s ... \n', mfilename);
printf('Using stimuli responses from %s\n', tools.file('T',stim_folder))
printf('Using membrane currents from %s\n', tools.file('T', AX.folder))

for f = fieldnames(EM.utils)' % Create utility local functions
  EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');
  eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
end

nerve = AX.nerve;

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

% Translate sensitivity peaks if specified
if any(named('-recenter-peak'))
  disp('Translating sensitivity peaks...')
  sensitivity = translate_I2V_peaks(sensitivity, get_('-recenter-peak'));
end

clear ee ff x_ y_ z_ fac_ f sel ok list

%% Load axon population, mem current index file and re-assign axon groups

list = dir([stim_folder filesep '*.mat']); 
load(p_(list(1)),'stimulus');

% idx = cellfun(@(n) str2double(regexp(n,'(?<=fascicle)\d+','match', ... 
%                                        'once')), {list.name},'unif',0);

has_fraction = ~isempty(whos('-file',p_(list(1)),'spike_fraction')); 

if has_fraction && ~any(named('-no-frac'))
    raster = arrayfun(@(n) load(p_(n), 'results','pop','axon_index','spike_fraction'), list);      
    fprintf('Using supplied downsampling (%0.1f-%0.1f%%)\n', ...
             100*min([raster.spike_fraction]), 100*max([raster.spike_fraction]))
else
    raster = arrayfun(@(n) load(p_(n), 'results','axon_index'), list);  
   [raster.spike_fraction] = deal(1); 
end

pop = AX.pop; 
[raster.filename] = deal(list.name);
full_raster = convert_raster_population(pop,raster); 



%% Get Distortion if specified

if any(named('-dis')),
  options.distortion = get_('-dis'); 
  if size(options.distortion,1) == numel(pop)
    options.distortion_by_type = options.distortion; 
  end
  % options.distortion. [ v z t ] is read as a matrix or structure by
  % models.spike_to_wave.
  
  if any(named('-debug-dis')), options.plot_image = true; end  
end

%%


fs = 30; % 24.414; % kHz 
if any(named('-sample-rate')), fs = get_('-sample-rate'); 
elseif any(named('-fs')),      fs = get_('-fs'); 
end

time_span = 65;
if any(named('-time')),  time_span = get_('-time'); end

n_stim = numel(stimulus.current); 
if ~any(size(stimulus.current) == 1)
    n_stim = size(stimulus.current,2); 
end
  
warn_once = true; 
options_setup = options; 

%% Compute intracellular-to-extracellular relationship

for i_stim = 1:n_stim

  waves = []; 
  raster = []; 
  options_list = []; 

  for ty = 1:numel(pop)
    time = 0:(1/fs):(time_span);
    time = [-fliplr(time) time(2:end)]; % ms

    options = options_setup; 
    options.filename = '';
    options.class = pop(ty).axon_model;
    options.padding = 2;
    options.spiketimes = [];

    % this is now handled in convert_raster_population
    xy = full_raster(ty).xy;
    grp = full_raster(ty).grp;
    fid = full_raster(ty).fid;
    gff = 10.^ceil(log10(nF+2)); % subgroup_id = gg + ff / 10 would fail for nF > 9

    if any(named('-fasc'))
      ok = ismember(fid,get_('-fasc'));
      xy = xy(ok,:); grp = grp(ok); fid = fid(ok); 
    end

    if isfield(options,'distortion_by_type')
      options.distortion = options.distortion_by_type(ty,:); 
    end

    % Work out whether the xy units need to be corrected
    if max(abs(xy(:))) > max(EM.model.nodes(:,1)),  xy = xy / 1e3;
        warning('ViNERS:axonUnits_mm',['%s should be in units of mm. ' ...
          'treating this as µm, maybe cancel this and fix?'], ...
          tools.file('T',[options.axons_folder pop(ty).axon_model '\index.mat']))

    elseif isfield(pop,'units') && isfield(pop.units,'axon_xy') && ...
          ~strncmp(pop.units.axon_xy,'mm',2)
      % the EIDORS data is in units of mm
      % the axon data /should/ also be in mm but has been known to pop
      % up in µm, which I think was an issue with fascicle tracings being
      % specified in µm.

      xy = xy / 1e3;
      warning('ViNERS:axonUnits_mm',['%s should be in units of mm, ' ...
        'not %s. Treating this as µm, maybe cancel this and fix?'], ...
        tools.file('T',[options.axons_folder pop(ty).axon_model{1} '\index.mat']),...
                     pop.units.axon_xy)

    elseif any(named('-unit-um')) xy = xy / 1e3;
    end % um -> mm ?

    nG = max(pop(ty).size_sample);

    options.efferent = false; % automatically determined from the input raster in this mode 

    S = full_raster(ty);
    S.dt_dx = [S.results.dt_dx]';
    S.spike = [S.results.spikes]';
    for aa = 1:numel(S.results)
        S.spike(aa).time = S.spike(aa).time{i_stim};
        S.spike(aa).node = S.spike(aa).node{i_stim};
        S.spike(aa).init = S.spike(aa).init{i_stim};
    end
    options.raster = rmfield(S,{'results'});
    options.raster.axon_group = grp + fid/gff; 
    
    if isempty(raster), raster = options.raster;
    else raster(ty) = options.raster; %#ok<AGROW>
    end

    if isempty(options_list), options_list = options;
    else options_list(ty) = options; %#ok<AGROW>
    end
    
    assert(numel(S.spike) == numel(grp),'population size mismatch')

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

    V = cell(nG,nF); % Observed potentials 

    if any(named('-debug-u')) || any(named('-preview'))
      %% Debug this, is a consistent problem ... 
      tools.setupEIDORS ;
      clf, show_fem(EM.model);
      hold on, h = get(gca,'Children');
      set(h,'EdgeAlpha',0.2); delete(h(1:2))

      for gg = 1:nG
        for uu = 1:size(axon_xy{gg},1)
          xyz = permute(axon_xy{gg}(uu,:,:), [2 3 1]);
          if size(axon_xy{gg},2) == 2
            plot3(0*xyz(1,:), xyz(2,:),xyz(1,:),'o', ...
                'Color',[1 0 0 0.2],'Clipping','off') 
          else
            plot3(xyz(:,1), xyz(:,2),xyz(:,3), ...
                'Color',[1 0 0 0.2],'Clipping','off') 
          end
          hold on
        end
      end
      
      view([1 0 0])
      %%
      return
    end

    printf('Computing spatial summation [%s], stimulus # %d ... ',pop(ty).axon_type,i_stim)

    %% Compute summation of currents from spike raster
    cache_path = tools.cache('path');

    if nG*nF == 1 || ~options.use_parallel
      for ii = 1:(nG*nF)
        V{ii} = parfun_unpack(cache_path, @(g,f) ...
                      models.spike_to_wave(g+f/gff,time,        ...
                                           sensitivity(:,ff),  ...
                                           axon_xy{g,f},options), ...
                                           [nG nF],ii);
      end
    elseif isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP> % octave parallel
        if isempty(which('pararrayfun')), pkg load parallel; end 
        V = pararrayfun(nproc-1, @(a) parfun_unpack(cache_path, @(g,f) ...
                                         models.spike_to_wave(g+f/gff,time, ...
                                                        sensitivity(:,ff), ...
                                                       axon_xy{g,f},options), ...
                                    [nG nF],a), 1:(nG*nF), 'Unif', false);
        error TODO_validate_this_for_octave
      else    
        parfor ii = 1:(nG*nF)  
          V{ii} = parfun_unpack(cache_path, @(g,f) ...
                        models.spike_to_wave(g+f/gff,time,        ...
                                             sensitivity(:,ff),  ...
                                             axon_xy{g,f},options), ...
                                          [nG nF],ii); %#ok<PFBNS>
        end
    end    
    printf('Done!\n')

    for ff = 1:nF % compute fascicle sum and apply spike-fration scaling
      waves(:,:,ff,ty) = sum(cat(3,V{:,ff}),3); %#ok<AGROW>
      waves(:,:,ff,ty) = waves(:,:,ff,ty) ./ S.spike_fraction(ff).^(sqrt(1/2)); 
    end
      
  end % ty [1..4]

  % This little check prevents from having to open another instance
  summary = [sqrt(nanmean(waves(:).^2)) nanmin(waves(:)) nanmax(waves(:))];
  printf(' computed wave: RMS %0.1f, range %f-%f\n', summary ); 
  if warn_once && all(summary == 0) && any([options_list.spiketimes])
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
      ax_g = (options.raster.axon_group == g_id); % Group mask
      ax_k = ax_g(options.raster.spk_axon); % which times are in group?
      if any(ax_k), ax_k = find(ax_k,1);      
        spk1_t(ii) = options.raster.spk_time(ax_k); 
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

    options.plot_wave = 1; 
    [gg,ff] = ind2sub([nG nF],ii); 
    [V_check,vSpk_check] = models.spike_to_wave(gg+ff/10, ... 
                                  time,sensitivity(:,ff),axon_xy{ii},options);

    %% Check the saved "results" field
    nG = numel(sam.A_diam);
    for gg = 1:nG % Load cache results 
      d = load(axonpath(ty,'n%03d_out.mat',gg),'result');      
      results(gg,:) = d.result;
    end

  end % if plot_debug

  %%
  if any(named('-f-sum'))   
    waves = squeeze(sum(waves,3)); 
  elseif size(waves,3) == 1, waves = squeeze(waves); 
  end

  clear gg ff ee g_xy ii xy grp ty fid ans

  if any(named('-out')), e_name = get_('-out');
  else
    % [~,e_name,~] = fileparts(stim_folder); % this mangles decimal numbers in the path
    e_name = regexp(stim_folder,'(?<=[\\/])[^\\/]+$','match','once');
    e_name = regexprep(e_name,'stimulus[^\w]*','');
    e_name = strtrim(regexprep(e_name,'\([^\)]*\)',''));
  end

  file_out = sprintf('stim_%03d.mat',i_stim);
  if isfield(stimulus,'CL')
      file_out = sprintf('stim_CL%03d.mat',stimulus.CL(i_stim));
  end
  wave_path = strrep(['stimulus (' e_name ')'],' ()','');
  wave_path = strrep(wave_path, '))',')');

  file_out = fullfile(tools.file('waves~/'),wave_path, ...
                      strrep(file_out,'.0_','_'));

  if ~exist(fileparts(file_out),'dir'), mkdir(fileparts(file_out)), end  

  options.class = {options_list.class};
  options.afferent = [pop.afferent];
  options.spiketimes = {options_list.spiketimes};
  options = rmfield(options,{'verbose','raster','efferent'});
  options.filename = file_out;  
  inputs = varargin; 
  
  printf('Saving %s\n', tools.file('TT',file_out))
  save(file_out, 'raster','waves','time','stimulus','inputs','options')

end % do_stimuli


if nargout == 0, clear, end
return

%% Move sensitivity peak (for electrode array property explorations) 
function sensitivity = translate_I2V_peaks(sensitivity, z_REF)
%%

if numel(z_REF) < nE, 
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

function R_full = convert_raster_population(pop,D)

% pop is the source axon population. D is a (nF*nTy) x 1 list with the
% spike-time, location data for each fascicle and axon model type. 

R = rmfield(D([]),'pop'); % initialise struct fieldnames
[R.population_id] = deal([]); 
types = unique(regexp({D.filename},'^[^-]*','match','once')); 

for ty = 1:numel(types)
  
  R(ty).filename = types{ty};  
  R(ty).results = D(1).results([]);  
    
  sel = find(contains({D.filename},types{ty}));   
  % assert(numel(sel) == nF,'Missing raster files')
  ok = false(max(cat(1,D(sel).axon_index)),1);
  
  for ii = sel  
      
    daxi = D(ii).axon_index;
    if any(ok(daxi))
      warning('ViNERS:duplicateSim', ...
              '%s Axon #%d (and %d others) have duplicate indices', ...
               types{ty}, find(ok(daxi),1), sum(ok(daxi))-1 )
    end
    
    ok(daxi) = true; 
    R(ty).results(daxi) = D(ii).results; 
    R(ty).axon_index(daxi) = D(ii).pop.source_index(daxi);
    R(ty).population_id(daxi) = D(ii).pop.population_id(daxi);
  end  
  
  R(ty).spike_fraction = [D(sel).spike_fraction];
  for f = {'results','axon_index','population_id'}      
    R(ty).(f{1}) = reshape(R(ty).(f{1})(ok),[],1);
  end  
end

% R is now a (1 x nTy) list with like models joined together across
% fascicles. 

%% break back apart into afferent/efferent populations  

D = R; 

for ty = 1:numel(pop)
  sel = contains({D.filename},pop(ty).axon_model);   
  R(ty) = D(sel);
  
  if length(unique(R(ty).population_id)) == 1, continue, end

  ok = (R(ty).population_id == ty);
  
  for f = {'results','axon_index','population_id'}
    R(ty).(f{1}) = R(ty).(f{1})(ok,:);
  end
  
  nF = size(R(ty).spike_fraction,2);
  R(ty).spike_fraction = arrayfun(@(ff) ... % recompute this
                       sum(pop(ty).fascicle(R(ty).axon_index) == ff) ./ ...
                       sum(pop(ty).fascicle == ff), 1:nF);  
end

[R.xy]  = deal([]); 
[R.grp] = deal([]);
[R.fid] = deal([]); 

for ty = 1:numel(pop)

    R(ty).xy = pop(ty).axon_xy(R(ty).axon_index,:);
    R(ty).grp = pop(ty).size_sample(R(ty).axon_index,:);
    R(ty).fid = pop(ty).fascicle(R(ty).axon_index,:);
end

R_full = R;

return


function wave = parfun_unpack(cache_path,fun,nGnF,ii)

  tools.cache('set',cache_path)
  [gg,ff] = ind2sub(nGnF,ii);
  wave = fun(gg,ff); 
