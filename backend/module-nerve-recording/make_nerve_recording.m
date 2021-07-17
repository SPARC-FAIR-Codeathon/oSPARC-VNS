
function make_nerve_recording( axons_file, fields_file, ...
                               spikes_file_or_string, sample_rate, ...
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

if nargin < 3 || isempty(spikes_file_or_string)    
  fprintf('Arg 2 not set, using default spikes settings: ')
  spikes_file_or_string = 'flat[0.2,2]';
  fprintf('%s\n', spikes_file_or_string)
end

if nargin < 4 || isempty(sample_rate)    
  sample_rate = 30; % kHz
  fprintf('Arg 4 not set, using default sample rate: %g ks/s\n', sample_rate)
end

named = @(v) strncmpi(v,varargin,length(v)); 
hasext = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 
tools.file('root',pwd); % set 'root' to this folder

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end

if true
    disp('====================================================')
    fprintf('Running models.nerve_recording %s\n', datestr(now))
    fprintf('{1} = %s\n{2} = %s\n{3} = %s\n',axons_file, fields_file, spikes_file_or_string)
    disp('====================================================')
end

if exist(tools.file('out~\waves'),'dir')
     rmdir(tools.file('out~\waves'),'s');
end, mkdir(tools.file('out~\waves')); 

inputs = [convert_eidors_file(fields_file) ... 
          unpack_axons_file(axons_file)]; 
      
%% Determine spiketimes to stimulate

if exist(spikes_file_or_string,'file')
  if hasext(spikes_file_or_string,'.mat')
    raster = load(spikes_file_or_string);
  elseif hasext(spikes_file_or_string,'.xml')
    raster = convert_xml_raster(spikes_file_or_string, sample_rate); 
  elseif hasext(spikes_file_or_string,'.json')
    raster = tools.parse_json(spikes_file_or_string);
  elseif hasext(spikes_file_or_string,'.nwb')
    raster = read_NWB_raster(spikes_file_or_string);
  else error('Unknown file format "%s", expected {.mat,.json,.xml}');
  end
  
else
  %% Parse spikes_file_or_string
  
  entries = regexp(spikes_file_or_string,'(flat|burst|drift)[^\]]+[^\,;\s]+','match'); 
  
  if isempty(entries)
    warning('no flat[], burst[], or drift[] entries provided to %s', mfilename)
    entries = {'flat[0.2,2]'};
  end
  
  coherence = 1; % default coherence
  repeats = 3;   % default # reps
  
  for ii = 1:numel(entries)
      
    nmr_setting = models.nerve_recording('-list',entries{ii}); 
    my_c = str2double(regexp(entries{ii},'(?<=_c)\d+(\.\d+)','match'));
    if ~isnan(my_c) && ~isempty(my_c), coherence = my_c; end
    my_r = str2double(regexp(entries{ii},'(?<=_c)\d+(\.\d+)'));
    if ~isnan(my_r) && ~isempty(my_r), repeats = my_r; end
    
    nmr_setting.coherence = coherence; 
    nmr_setting.n_reps = repeats;
    
    input_args = regexp(entries{ii},'(?<=\[)[^\]]+','match','once');
    
    switch(lower(entries{ii}(1:4)))
     case 'flat'
       nmr_setting.spikerate = reshape(str2double(strsplit( ...
                                       input_args,{',',';',' '})), [], 1); 
     case 'burs'
         
         
         
         
     case 'drif'
         
         
         
    end
  end
  
end





models.nerve_recording(inputs{:})

disp('TODO: condense epochs into a single file to return to user')

return







%% generate expected file structure in tempdir
function args = unpack_axons_file(filename) 


load(filename,'axons','axon_models','nerve');

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
save(filename,'-struct','EM')
args = {'-file', filename};
  
return

function raster = convert_xml_raster(filename, fs)


%%
if ~exist('filename','var'), filename = 'test-raster.xml'; end

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
info.Data.filename = filename;
info.Data.syntax = xml.Name;

for ii = 1:numel(info.Children)
  val = attr_(info.Children(ii),'value'); 
  if ~isnan(str2double(val)), val = str2double(val); end
  info.Data.(attr_(info.Children(ii),'name')) = val; 
end, info = info.Data; 

%% Get spike-times

spk_time = []; 
spk_unit = []; 

spikes = the(xml,'s');

for ii = 1:numel(spikes)    
  spk_time = [spk_time; str2double(attr_(spikes(ii),'t'))]; %#ok<AGROW>
  spk_unit = [spk_unit; str2double(attr_(spikes(ii),'u'))]; %#ok<AGROW>
end

units = the(xml,'unit');

for ii = 1:numel(units)

    uid = str2double(attr_(units(ii),'id'));    
    t = textscan(units(ii).Children(1).Data,'%f');

    spk_time = [spk_time; t{1}];       %#ok<AGROW>
    spk_unit = [spk_unit; 0*t{1}+uid]; %#ok<AGROW>
end

ok = ~isnan(spk_time) & ~isnan(spk_unit); 
if ~all(ok)
    
  warning('In %s, %d/%d spiketimes were corrupted.', ...
            filename, sum(~ok), numel(ok))
  spk_time = spk_time(ok); 
  spk_unit = spk_unit(ok);
end


if isfield(info,'time_start') && isfield(info,'time_end')  
     time = (info.time_start) : (1/fs) : (info.time_end);   
else time = floor(min(spk_time)) : (1/fs) : ceil(max(spk_time)); 
end

if isfield(info,'n_units'), nA = info.n_units;
elseif isfield(info,'n_axons'), nA = info.n_axons;
elseif isfield(info,'n_cells'), nA = info.n_cells;
elseif isfield(info,'n_spikes'), nA = info.n_spikes;
else nA = max(spk_unit); 
end

if isfield(info,'bin_dt')    
     bin_x = min(time):(info.bin_dt):max(time);
else bin_x = min(time):(1):max(time); % default 1 ms bins
end, bin_y = hist(spk_time,bin_x);  %#ok<HIST>
bin_y = bin_y ./ mean(diff(bin_x)) / nA * 1000;

raster.spk_time = spk_time;
raster.spk_axon = spk_unit;
raster.spk_rate = bin_y;
raster.bin_time = bin_x; 
raster.bin_rate = []; % undefined
raster.pop_rate = numel(spk_time) / range(time) / nA * 1000;

raster.file_info = info; 

return


function generate_NWB_raster(inputs)




error here


nwb = NwbFile( ...
    'session_description', 'o2S2PARC simulation of nerve recording',...
    'identifier', 'Mouse5_Day3', ...
    'session_start_time', datetime(now,'ConvertFrom','datenum'), ...    
    'general_session_id', 'session_1234', ... % optional
        'ghjkl' )
    
    
trials = types.core.TimeIntervals( ...
    'colnames', {'start_time', 'stop_time', 'correct'}, ...
    'description', 'trial data and properties', ...
    'id', types.hdmf_common.ElementIdentifiers('data', 0:2), ...
    'start_time', types.hdmf_common.VectorData('data', [.1, 1.5, 2.5], ...
   	'description','start time of trial'), ...
    'stop_time', types.hdmf_common.VectorData('data', [1., 2., 3.], ...
   	'description','end of each trial'), ...
    'correct', types.hdmf_common.VectorData('data', [false, true, false], ...
   	'description', 'whether the trial was correct'))
    


    
function data = r
    error('TODO: understand how to read NWB files')
    
    
    
    
    
    
    
    
    