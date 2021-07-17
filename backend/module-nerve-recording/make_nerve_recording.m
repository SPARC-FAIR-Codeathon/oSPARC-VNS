
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
  spikes_file_or_string = 'flat[1,2,5,10]';
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


inputs = [convert_eidors_file(fields_file) ... 
          unpack_axons_file(axons_file)] 

if exist(spikes_file_or_string,'file')
  if hasext(spikes_file_or_string,'.mat')
    raster = load(spikes_file_or_string);
  elseif hasext(spikes_file_or_string,'.xml')
      
  elseif hasext(spikes_file_or_string,'.xml')
      raster = load(spikes_file_or_string);       
  else error('Unknown file format "%s", expected {.mat,.json,.xml}');
  end
    
else
    
  disp TODO_pre_save_NWB_rasters
end
% 
% 
% 
% if ~isempty(nerve_script)
%  
%   printf('Loading %s', nerve_script);
%   if has_ext_(nerve_script,'.mat'), n = load(nerve_script); 
%   elseif has_ext_(nerve_script,'.json'), n = tools.parse_json(nerve_script);
%   else error('%s: unknown extension, expected .mat or .json', nerve_script)
%   end      
% 
%   if isfield(n,'nerve'), opts.nerve = n.nerve; else opts.nerve = n; end
%   if isfield(n,'mesh'), opts.mesh = n.mesh; end
% end
% 
% if ~isempty(nerve_xml), opts.nerve.file = nerve_xml;
% elseif ~isfield(opts,'nerve'), opts.nerve.file = nerve_xml;
% end


% mnr_args = {'-coh',c,'-reps',r}


models.nerve_recording(inputs{:})

error TODO







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
    if sum(sel) ~= 1,
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

