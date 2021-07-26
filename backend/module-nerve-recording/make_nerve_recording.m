
function make_nerve_recording( axons_file, fields_file, ...
                               spikes_string, spikes_file, sample_rate, ...
                               varargin)
% Make_nerve_recording

if nargin < 1 || isempty(axons_file)
  fprintf('Arg 1 not set, using demo axons.mat file \n')
  axons_file = './input/demo/axon-population (1).mat'; 
end

if nargin < 2 || isempty(fields_file)
  fprintf('Arg 2 not set, using demo file extracellular-potential.mat file \n')
  fields_file = './input/demo/extracellular-potential (1).mat'; 
end

if nargin < 3, spikes_string = 'flat[0.2,2]';
  fprintf('Arg 3 not set, using default spikes settings: ')
  fprintf('%s\n', spikes_string)
end

if nargin < 4 || (isempty(spikes_string) && isempty(spikes_file))
    spikes_file = ''; 
  % spikes_file = './input/demo/example.xml'; 
  % fprintf('Arg 4 not set, using demo raster xml file: ')
  % fprintf('%s\n', spikes_file)
elseif ~isempty(spikes_string) && ~contains(spikes_string,'file')
  spikes_file = ''; 
end

if nargin < 5 || isempty(sample_rate)    
  sample_rate = 30; % kHz
  fprintf('Arg 5 not set, using default sample rate: %g ks/s\n', sample_rate)
elseif ischar(sample_rate), sample_rate = str2double(sample_rate); 
end

if isdeployed, warning('off','BIDSfile:notFound'), end

named = @(v) strncmpi(v,varargin,length(v)); 
tools.file('root',pwd); % set 'root' to this folder


if true
    disp('====================================================')
    fprintf('Running models.nerve_recording %s\n', datestr(now))
    fprintf('{1} = %s\n{2} = %s\n{3} = %s%s\n',axons_file, fields_file,  ...
                spikes_string, spikes_file)
    disp('====================================================')
end

if exist(tools.file('out~\waves'),'dir')
     rmdir(tools.file('out~\waves'),'s');
end, mkdir(tools.file('out~\waves')); 

inputs = [convert_eidors_file(fields_file) ... 
          unpack_axons_file(axons_file), '-fs', sample_rate]; 
      
if isdeployed, disp('Progress: 5%'); end
      
if any(named('-:')), inputs = [inputs varargin(find(named('-:'))+1:end)]; 
end


%% Determine spiketimes to stimulate

if exist(spikes_file,'file')
  models.nerve_recording(inputs{:}, '-raster', spikes_file )
else
  %% Parse string command
  
  entries = regexp(spikes_string,'(flat|burst)[^\]]+[^\,;\s]+','match'); 
  
  if isempty(entries)
    warning('no flat[] or burst[] entries provided to %s', mfilename)
    entries = {'flat[0.2,2]'};
  end
  
  coherence = 1; % default coherence
  repeats = 3;   % default # reps
  
  for ii = 1:numel(entries)
      
    mnr_setting = models.nerve_recording('-list',entries{ii}); 
    my_c = str2double(regexp(entries{ii},'(?<=_c)\d+(\.\d+)?','match'));
    if ~isempty(my_c) && ~any(isnan(my_c)), coherence = my_c; end
    my_r = str2double(regexp(entries{ii},'(?<=_r)\d+','match'));
    if ~isempty(my_r) && ~any(isnan(my_r)), repeats = my_r(1); end
    
    mnr_setting.coherence = coherence; 
    mnr_setting.n_reps = repeats;
    
    input_args = regexp(entries{ii},'(?<=\[)[^\]]+','match','once');
    input_args = str2double(strsplit(input_args,{',',';',' '})); 
    switch(lower(entries{ii}(1:4)))
     case 'flat'
       mnr_setting.spikerate = reshape(input_args, [], 1); 
     case 'burs'
               
       mnr_setting.name = 'gauss';          
       mnr_setting.spikerate = [input_args(1:3:end); input_args(2:3:end)]';
       mnr_setting.frequency =  input_args(3:3:end)/sqrt(8*log(2)); 
       mnr_setting.exponent = 1; 
       mnr_setting.wave_path = 'burst';
       mnr_setting.file_scheme = 'epoch_k%0.1f_w%0.1f_c%0.1f (%%d).mat';
       mnr_setting.file_vector = @(ex,sr,fr,ch) [sr(2) fr ch];
    end
    
    fprintf('Set intra-population coherence = %g\n', coherence);
    fprintf('Set experimental replicates = %g\n', repeats);
    
    if isdeployed, fprintf('Progress: %0.1f%%\n', 5+90*(ii-1)/numel(entries)); end
    
    mnr_setting.wave_path = sprintf('%s_%02d', mnr_setting.wave_path, ii);
    models.nerve_recording(inputs{:}, '-settings', mnr_setting )

  end
  
end

if isdeployed, disp('Progress: 95%'); end

merge_stimulation_epochs(inputs)

if isdeployed, disp('Progress: 100%'); end
return







%% generate expected file structure in tempdir
function args = unpack_axons_file(filename) 

fprintf('Loading %s ...\n', filename)
load(filename,'axons','axon_models','nerve');

if isempty(axons(end).fascicle)
  disp('warning: buggy axons(4).fascicles')
  n4 = numel(axons(4).fibre_diam)-1;
  axons(4).axon_xy = axons(3).axon_xy(end-n4:end,:);
  axons(4).fascicle = axons(3).fascicle(end-n4:end);
  axons(3).axon_xy(end-n4:end,:) = []; 
  axons(3).fascicle(end-n4:end) = []; 
end

model_names = unique({axons.axon_model});

if exist(tools.file('in~\axons'),'dir')
     rmdir(tools.file('in~\axons'),'s');
end, mkdir(tools.file('in~\axons')); 

% cache = @(varargin) tools.cache('path',sprintf(varargin{:}));
cache = @(varargin) ['./input/' sprintf(varargin{:})];

for type = model_names
  
  mkdir(cache('axons/in/%s',type{1}))
    
  for gg = 1:numel(axon_models.(type{1}).simulation)
    this = axon_models.(type{1}).simulation(gg);
    save(cache('axons/in/%s/n%03d_out.mat',type{1},gg), '-struct', 'this')
  end  
  
  [pop,sam] = merge_populations(axons,type);
  results = axon_models.(type{1}).summary; 
  save(cache('axons/in/%s/index.mat',type{1}), 'pop','sam','results')
end

pop = axons; 
opts = {'generated using unpack_axons_file'}; 
save(cache('axons/axons (in).mat'),'pop','opts','nerve');

args = {'-axon-file',cache('axons/axons (in).mat'), ...
        '-currents',cache('axons/in/')}; 

return



%% All axon populations with a certain axon model
function [pop,sam] = merge_populations(pop,axon_type)

  persistent sam_myel sam_unmy
  if ischar(pop), sam_myel = []; sam_unmy = []; return, end
  nG = max(cat(1,pop.size_sample)); % # of sub-groups

  if isempty(sam_myel)      
    %% construct sample across all myelinated axons
    sam_myel = struct;
    sel = [pop.myelinated]; 

    sam_myel.axon = cat(1,pop(sel).size_sample);
    sam_myel.count = arrayfun(@(x) sum(sam_myel.axon == x), 1:nG)';

    the = @(v) cat(1,pop(sel).(v)); 
    sam_myel.fibre_diam = the('fibre_diam');
    sam_myel.axon_diam = the('axon_diam');
    sam_myel.g_ratio = the('g_ratio');
    
    the = @(v) arrayfun(@(x) median( v(sam_myel.axon == x)), (1:nG)');
    sam_myel.fibre_diam = the(sam_myel.fibre_diam);
    sam_myel.axon_diam = the(sam_myel.axon_diam);
    sam_myel.g_ratio = the(sam_myel.g_ratio);

    sam_myel.subtype_id = find(sam_myel.count);
    
    if any(sam_myel.count == 0)
      warning('ViNERS:membraneCurrent:emptySubgroups', ...
              'The sub-sampling for %s had %d / %d empty sample groups.', ...
              'myelinated axons', sum(sam_myel.count == 0), numel(sam_myel.count))
      
      missing = (sam_myel.count == 0);
      for f = fieldnames(sam_myel)'
        if numel(f) ~= numel(missing), continue, end
        sam_myel.(f{1})(missing) = []; 
      end
    end
    
    %% construct sample across all unmyelinated axons axons
    
    sam_unmy = struct;
    sel = ~[pop.myelinated]; 

    sam_unmy.axon = cat(1,pop(sel).size_sample);
    sam_unmy.count = arrayfun(@(x) sum(sam_unmy.axon == x), 1:nG)';

    the = @(v) cat(1,pop(sel).(v)); 
    sam_unmy.fibre_diam = the('fibre_diam');
    
    the = @(v) arrayfun(@(x) median( v(sam_unmy.axon == x)), (1:nG)');
    sam_unmy.fibre_diam = the(sam_unmy.fibre_diam);
    
    sam_unmy.subtype_id = find(sam_unmy.count);

    if any(sam_unmy.count == 0)
      warning('ViNERS:membraneCurrent:emptySubgroups', ...
              'The sub-sampling for %s had %d / %d empty sample groups.', ...
              'unmyelinated axons', sum(sam_unmy.count == 0), numel(sam_unmy.count))
      
      missing = (sam_unmy.count == 0);
      for f = fieldnames(sam_unmy)'
        if numel(sam_unmy.(f{1})) ~= numel(missing), continue, end
        sam_unmy.(f{1})(missing) = []; 
      end
    end
  end

  [pop.population_id] = deal([]);
  for ty = 1:numel(pop)
      pop(ty).population_id = ty * ones(size(pop(ty).size_sample)); 
  end
  sel = strcmp({pop.axon_model},axon_type{1}); 
  pop = pop(sel);   
  
  for f = fieldnames(pop)'
    if numel(pop) <= 1, break, end
    if ischar(pop(1).(f{1})), pop(1).(f{1}) = {pop.(f{1})};
    else pop(1).(f{1}) = cat(1,pop.(f{1})); end
  end
  pop = pop(1); pop(1).axon_model = axon_type{1};
  assert(all(pop.myelinated) == any(pop.myelinated),'Myelination is all-or-none')
  pop.myelinated = pop.myelinated(1);
    
  if pop.myelinated, sam = sam_myel; else sam = sam_unmy; end

%% Convert v_extracellular to structure 
function args = convert_eidors_file(filename)

fprintf('Loading %s ...\n', filename)
EM = load(filename); 
% see https://en.wikipedia.org/wiki/Reciprocity_(electromagnetism)#Reciprocity_for_electrical_networks               

nF = sum(contains(EM.model.object_name,'Fasc')) - ...
     sum(contains(EM.model.object_name,'P_Fasc'));

fasc_ = @(n) sprintf('Fascicle%d',n);

% EM.v_extracellular is a list, (nodes x electrodes) of the sequential
% bipolar stimulus-induced extracellular potential, for ALL nodes 
% (not just fascicle nodes). 

% a sensitivity file is broken up by fascicle, which each fascicle
% containing a NODE index and a list of per-node values 

% EM.model.object_id is a list of ELEMENT ids which make up each object
% in the EIDORS model. Gather each NODE in that set of elements, then use
% that to emulate the EM.FascicleN.pot/idx 

for ff = 1:nF

    sel = strcmpi(EM.model.object_name,fasc_(ff));
    if sum(sel) == 0, sel = strcmpi(EM.model.object_name,fasc_(1)); end
    if sum(sel) ~= 1
        error('Object %s not found in EM.model.object_name',fasc_(ff)); 
    end
    idx = unique(EM.model.elems(EM.model.object_id{sel},:));

    EM.(fasc_(ff)).idx = idx'; 
    EM.(fasc_(ff)).pot = EM.v_extracellular(idx,:)';
end  

EM.utils.x_ = '@(i) out.model.nodes(i,1)';
EM.utils.y_ = '@(i) out.model.nodes(i,2)';
EM.utils.z_ = '@(i) out.model.nodes(i,3)';
EM.utils.fac_ = '@(i) out.(sprintf(''Fascicle%d'',i)).idx'; 

EM = rmfield(EM,'v_extracellular');
EM.info.filename = filename; 

if exist(tools.file('in~\eidors'),'dir')
     rmdir(tools.file('in~\eidors'),'s');
end, mkdir(tools.file('in~\eidors')); 


filename = tools.file('in~\eidors\sensitivity.mat'); 
fprintf('Saving %s ...\n', filename)

save(filename,'-struct','EM')
args = {'-file', filename};
  
return


function merge_stimulation_epochs(inputs)

list = dir('./output/waves/**/*.mat'); 
p_ = @(x) [x.folder filesep x.name]; % path expander

load(inputs{4},'pop','nerve');

data = struct;
data.anatomy = nerve; 
data.epoch_labels = {};
data.population = []; 
data.options = []; 
data.time = []; 
data.input_arguments = []; 

for ii = 1:numel(list)
  
  d = load(p_(list(ii)));
  data.epoch_labels{ii,1} = list(ii).name(1:end-4);
  
  roi = (abs(d.time) <= max(d.time)-15);  
  d.waves = tools.detrend_wave(d.waves,d.time,roi); 
  d.waves = d.waves(roi,:,:,:);
  d.time = d.time(roi);
  
  for ty = 1:numel(d.raster)
     
    if numel(data.population) < ty,
       
      this = struct; 
      this.class = pop(ty).axon_type;
      this.axon_info = pop(ty);      
      
      this.waves{1} = permute(d.waves(:,:,ty,:),[1 2 4 3]); 
      this.spikes    = d.raster(ty);
      
      if isempty(data.population), data.population = this;
      else data.population(ty,1) = this;
      end
    else
      data.population(ty).waves{ii,1} = permute(d.waves(:,:,ty,:),[1 2 4 3]);
      data.population(ty).spikes(ii,1) = d.raster(ty);
    end
  end
  
  if isempty(data.options)
       data.options = d.options; 
       data.time = d.time;
       data.input_arguments = d.inputs;
  else data.options = [data.options; d.options];
  end
end
    

for ty = 1:numel(d.raster)
    
  data.population(ty).waves = cat(4,data.population(ty).waves{:});
  
  data.population(ty).waves = permute(data.population(ty).waves, [4 1 2 3]);
end
    

data.units.waves = 'micro volt';
data.units.time = 'milli second';
data.units.spikerate = 'impulse per second per axon';
data.units.waves_dimensions = {'epoch','time','electrode','fascicle'};
data.units.waves_configuration = 'monopolar recording';

filename = tools.file('get','out~/nerve-recording (%d).mat','next');

fprintf('Combining waves files into %s\n', filename)
save(filename,'-struct','data','-v7.3')

    
    
    
    
    
    
    
    