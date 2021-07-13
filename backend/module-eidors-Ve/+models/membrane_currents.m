
function membrane_currents(varargin)
% models.membrane_currents computes the membrane currents and conduction
% velocites for an axon population ( specified in axons~/axons.mat ). 
% 
% The following arguments are supported: 
%  -file     : specify input axon populations (default: newest axons.mat)
%  -out      : specify output folders location
%  -no-pct   : run in debug mode (no parallel compute). 
%  -args {}  : extra arguments to models.axon_model
%
%  -fix-Gaines, -fix-Sundt, -fix-MRG : re-run just one class 
% 
% v0.3 CDE 13 June 2021
% v0.2 CDE 30 May 2020

varargin = tools.opts_to_args(varargin,'axons');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

workingDir = pwd; 
cd(tools.file('~/code/'))

% get source filename
if any(named('-file')), index_file = get_('-file'); 
else index_file = tools.file('get','axons~/axon*.mat','newest'); 
end
if any(index_file == '~'), index_file = tools.file(index_file); end
if ~exist(index_file,'file')
    index_file = tools.file('get',['axons~/*' index_file '*.mat']);    
end

disp(['Loading ' tools.file('short',index_file)])
load(index_file,'pop');

FULL_axon_population = pop; 
merge_populations('reset'); 

axon_type_list = unique({pop.axon_model});
output_stub = regexp(index_file,'(?<=\()[^\)]+','match','once'); 

for axon_type = axon_type_list
  
  if any(named('-fix-'))
    if ~any(named(['-fix-' axon_type{1}])), continue, end
  end
  if ~any(named('-debug-review')) % Can skip if needed for debugging
    tools.cache('reset');
  end
  
  %% Extract relevent population 
  
  [pop,sam] = merge_populations(FULL_axon_population,axon_type);
  
  nG = numel(sam.subtype_id); % # of sub-groups 

  
  %% Work out destination folder for this axon_type
  if any(named('-out')), out_folder = get_('-out'); 
    if any(out_folder == '~'), out_folder = tools.file(out_folder);
    elseif ~any(ismember(out_folder,'/\')) 
      out_folder = [tools.file('sub~/axons/') out_folder '/']; %#ok<AGROW>
    end
  else out_folder = tools.file(['sub~/axons/',output_stub,'/']);
  end, out_folder = sprintf('%s%s%s/', out_folder, axon_type{1}); 
  
  if exist(out_folder,'dir') 
    copyfile(out_folder,tools.file('get', [out_folder ...
                          '/../duds/' axon_type{1} '(%d)'],'next'),'f')
    rmdir(out_folder,'s'); 
  end, mkdir(out_folder)
  out_folder = regexprep(out_folder,'[\\/]+',filesep); 
  save([out_folder '/index.mat'],'pop');

  %% Run a membrane current output for each sample of axons
  % nG = length(sam.count); 

  mdl_args = {[],pop.axon_model,'-current','-constantLength',8};
  if any(named('-debug')), mdl_args = [mdl_args {'-debug'}]; end           %#ok<AGROW>
  if any(named('-arg')), extra_args = get_('-arg'); 
    if ~iscell(extra_args), model_args = [model_args extra_args];          %#ok<AGROW>
    else model_args = [model_args extra_args{:}];                          %#ok<AGROW>
    end
  end
  
  if pop.myelinated
      run_model = @(g) models.axon_model(sam.subtype_id(g),mdl_args{:}, ...
                               'fibreDiam',sam.fibre_diam(g), ...
                                 'g_ratio',sam.g_ratio(g));
  else mdl_args{5} = 6; % default length 6 mm
      run_model = @(g) models.axon_model(sam.subtype_id(g),mdl_args{:}, ...
                               'fibreDiam',sam.fibre_diam(g));      
  end

  fprintf('runing I_m NEURON model (%s) ... \n',axon_type{1})
  cache_path = tools.cache('reset');
  results = cell(nG,1); 

  %% CORE parallel loop to execute models.axon_model
  
  if any(named('-no-p'))
    for gg = 1:nG, results{gg} = feval(run_model,gg); end     
  elseif isempty(strfind(ctfroot, 'MATLAB')) %#ok<STREMP> % octave parallel
    if isempty(which('pararrayfun')), pkg load parallel; end
    
    results = pararrayfun(nproc-1, @(a) { tools.cache('set',cache_path), ...
                                          run_model(a)}, axons, ...
                                          'UniformOutput',0);
    results = [results{:}];  % comes out as 1x2 cell array
    results = [results{2:2:end}]; % convert to struct array      
  else                     
    set_cache = @() tools.cache('set',cache_path); 
    parfor gg = 1:nG,
      set_cache(); % otherwise files get lost
      results{gg} = feval(run_model,gg); 
    end
  end
  
  list = dir(tools.cache('path','n*.mat'));
  list(1).gid = str2double(regexp({list.name},'\d+','match','once'));
    
  for gg = sam.subtype_id' % file existance check
    if any(list(1).gid == gg), continue, end  
    fprintf('[%03d] Patching ... \n', gg)
    if ~ischar(run_model), run_model = func2str(run_model); end
    if strcmpi(axon_type{1},'Sundt'), hi_range = [2 1 10]; % from 0.25 0 2
    else                              hi_range = [1 0 2]; % from 0.05 0 0.5
    end
    eval(sprintf('p_model = %s,''-sweep-range'',[%g %g %g]);', ...
                    run_model(1:end-1), hi_range))
    results{gg} = feval(p_model,gg);
  end
  
  %% Collect results 
  
  results = [results{:}]; % convert to struct array
  if ~all(sam.subtype_id' == 1:nG) % if any empty groups  
     results(sam.subtype_id) = results;
     missing = setdiff(1:numel(results),sam.subtype_id);
    [results(missing).threshold] = deal(nan);
    [results(missing).spiketime] = deal(nan);
    [results(missing).velocity] = deal(nan);
  end
  
  list = dir(tools.cache('path','n*.mat')); % get any extras created during Patching...
  list(1).gid = str2double(regexp({list.name},'\d+','match','once'));
  
  for gg = 1:numel(list)
    % move file (and convert from 0000x to 00x)
    cmd = sprintf('move "%s%s%s" "%s/n%03d_out.mat"', ... 
              list(gg).folder,filesep,list(gg).name,out_folder,list(1).gid(gg));
    cmd = regexprep(cmd,'[\\/]',filesep); % "MOVE" is picky about slashes
    system(cmd);
  end
  
  fprintf('\n%s Membrane current profiles calculated.\n',axon_type{1})
  save([out_folder '/index.mat'],'pop','sam','results');

end

%% 

cd(workingDir)
return
%% Patch results 

disp('Patching ... ') %#ok<UNRCH>

fin = dir([out_folder filesep '*_out.mat']); 
fin = str2double(regexp({fin.name},'\d+','match','once')); 
missing = find(~ismember(1:size(results,1),fin)); 
% missing = find(isnan(results(:,3)));
% missing = find(results(:,2) < 0.001);

mresults = results(missing,:); 

parfor u = 1:numel(missing), ii = missing(u);
  mresults(u,:) = feval(run_model,ii); %%#ok<PFBNS>
end

for u = 1:numel(missing), ii = missing(u);
  system(sprintf('move "%sn%03d_out.mat" "%s"',tools.cache('path'),ii,out_folder));
  % models.axon_sfap and possibly others assume %03d label
end

results(missing,:) = mresults;
% results(isnan(results(:,1)),1) = inf; 

clear u fin missing mresults



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


  
  

