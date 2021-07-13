
function goal_metric = nerve_stimulation(varargin)
% models.nerve_stimulation(eidors_file, varargin) runs predefined stimulus
% levels for each for each axon in the sample (defined in sub~/axons/axons.mat). 
% 
% The following options are supported:
%  -pot (eidors~/stimulus.mat) : source file for extracellular stimulation
%           potentials (if not specified as the leading argument) 
%  -stim : which stimulus in eidors~/stimulus.mat to use? 
%  -axon (axons~/axons.mat) : specifiy axon file to simulate. 
%  -xy : run an XY test, ignoring axon population information. 
%  -resume : The simulations can take a long time. if -resume is set,
%         existing files will be skipped (e.g. you are resuming a sim)
% -fix-Gaines, -fix-Sundt, -fix-MRG : run just the specified class
%
% Adapted from models.threshold_model 
varargin = tools.opts_to_args(varargin,'stimulation');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

%% Parse input arguments 

working_dir = pwd; 
if any(named('-q')), printf = @(varargin) []; else printf = @fprintf; end
printf('Running models.%s ... \n', mfilename);

if isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP>,
    save_default_options ('-mat-binary')
end

[EM,AX] = tools.parse_arguments(varargin, 'LOAD', ...
                              'eidors','eidors~/stim*.mat', 'axons');

[full_population,input_stimulus] = process_input_args(get_,named,EM,AX);                         
                          
% Set the output path (Default: sub/stimulation/(eidors filename)
output_path = EM.filename; 
if any(named('-out')), output_path = get_('-out'); end
if isa(output_path,'function_handle'), output_path = output_path(); end
if ~any(ismember(output_path,'/\')), % if not already a path ... 
  output_path = tools.file(['stim~/' output_path '/'],'-make');
end
if any(named('-axon-t')), EM.info.AxonTrajectory = get_('-axon-t'); end

printf('Output will be saved to %s[...]\n', tools.file('T',output_path));
axon_type_list = unique({full_population.axon_model});
cd(tools.file('~/code/'))

%% MAIN loop

for ff = fascicle_list % Loop over fascicles 

  for axon_type = axon_type_list
          
    if any(named('-fix-')) % go back and fix one of them
      if ~any(named(['-fix-' axon_type{1}])), continue, end
    end
    
    %%
    stimulus = input_stimulus;
    pop = merge_populations(full_population,axon_type);     

    
    if isfield(stimulus,'name'), stimname = ['-' stimulus.name];
    elseif isfield(stimulus,'id'), stimname = sprintf('-stim%02d',stimulus.id);
    else stimname = '';
    end
    
    output_file = sprintf('%s%s-fascicle%d%s.mat', output_path, ... 
                                                 axon_type{1}, ff, stimname);
    if any(named('-xy')), 
      output_file = strrep(output_file,'.mat','-xy.mat');
    end

    if exist(output_file,'file') && any(named('-resume')), continue, end
    cache_path = tools.cache('reset');

    if any(named('-Ve')), Ve = get_('-Ve'); % this can be used to pass in other Ve patterns (e.g. constructued stimuli)
      if iscell(Ve), Ve = Ve{ff}; end
    else Ve = tools.load_Ve_field(EM,'-fascicle',ff, ... % Get Ve from EM 
                                     '-stim',1:size(stimulus.p,2),'-all'); 
    end
    
    model_args = {Ve,pop.axon_model,'-spikes','-constantLength',8, ...
                     '-fixed-stimulus-levels','-stimulus',stimulus};
    if any(named('-volt')), model_args{3} = '-voltage'; end
    if any(named('-debug')), model_args{end+1} = '-debug'; end %#ok<AGROW>
    if any(named('-arg')), extra_args = get_('-arg'); 
      if ~iscell(extra_args), extra_args = {extra_args}; end
      model_args = [model_args extra_args]; %#ok<AGROW>
    end

    if tools.from_trajectory(EM, AX.nerve), 
         axon_xy = tools.from_trajectory(EM, AX.nerve, pop.axon_xy);
    else axon_xy = pop.axon_xy;
    end
    
    if pop.myelinated
       run_model = @(g) models.axon_model(g,model_args{:}, ...
                                   'fibreDiam',pop.fibre_diam(g), ...
                                     'g_ratio',pop.g_ratio(g), ...
                                         '-xy',axon_xy(g,:,:));
    else model_args{5} = 6; % default length 6 mm
       run_model = @(g) models.axon_model(g,model_args{:}, ...
                                   'fibreDiam',pop.fibre_diam(g), ...
                                         '-xy',axon_xy(g,:,:));
    end
    
    axon_index = find(pop.fascicle == ff);    
    nA = numel(axon_index);

    if any(named('-repair-axon-indices'))
       fprintf('Repairing axon indices, %s Fascicle%d (%d axons)\n', ...
                          axon_type{1}, ff, nA)

        if ~exist(output_file,'file')
            warning('%s missing', output_file)
            continue
        end
                        
        d = load(output_file);
        d.pop = pop;
        d.axon_index = axon_index;
        save(output_file,'-struct','d')
        continue
    end
    
    fprintf('\n%s\nSimulating %s Fascicle%d (%d axons)\n', ...
                          '='*ones(1,40),axon_type{1}, ff, nA)
                      
    cache_path = tools.cache('reset');     
    subs_info = tools.file('info'); 
    set_cache = @() tools.cache('set',cache_path);
    results = cell(nA,1);
    
    %% CORE : parallel execution of models.axon_model

    
    
    if any(named('-no-p'))
      for aa = 1:nA, results{aa} = feval(run_model,axon_index(aa)); end
      results = [results{:}]; % comes out as Nx1 cell array
      
    elseif isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP> % octave parallel
      if isempty(which('pararrayfun')), pkg load parallel; end
      
      results = pararrayfun(nproc-1, @(a) { set_cache(), ...
                                            run_model(a)}, axon_index, ...
                                            'UniformOutput',0);
      results = [results{:}];  % comes out as 1x2 cell array
      results = [results{2:2:end}]; % convert to struct array      
    else       
      parfor aa = 1:nA,
        set_cache(); % otherwise files get lost
        results{aa} = feval(run_model,axon_index(aa)); 
      end
      results = [results{:}]; % comes out as Nx1 cell array      
    end

    %% Format output
    
    threshold = [results.threshold]'; %#ok<*NASGU> % (:,1); 
    velocity = [results.velocity]'; % results(:,2); 
    diameter = pop.fibre_diam(axon_index);
    
    if isfield(results,'spiketime'), spiketimes = {results.spiketime}'; end
    
    if any(named('-downs'))
      src_total = sum(pop.source_total(:,ff)); 
      spike_fraction = numel(threshold) / src_total;
      % axon_index = pop.source_index(axon_index); 
    end
    
    if ~isempty(pop.g_ratio), g_ratio = pop.g_ratio(axon_index); end
       
    clear results     
    list = dir(tools.cache('path','n*_out.mat'));
    [~,seq] = sort(str2double(regexp({list.name},'\d+','match','once'))); 
    list = list(seq);
    
    v_list = {'dt_dx','spikes','length','result'};
    results = arrayfun(@(a) load([a.folder filesep a.name],v_list{:}),list);
    n_nodes = arrayfun(@(r) length(r.length), results,'unif',0);
   [results.summary] = deal(results.result);
   [results.n_nodes] = deal(n_nodes{:});
   
    spike_counts = cellfun(@(s) cellfun(@(n) size(n,1), s.init), {results.spikes},'unif',0);
    spike_counts = cat(1,spike_counts{:});
   
    results = rmfield(results,v_list(4:end));
                          
    units.time = 'ms';
    units.length = 'mm';
    units.diameter = 'um';
    units.threshold = 'uA';
    units.stim_pattern = 'relative';
    units.stim_current = 'uA';
    units.spikes = {'axon_id','stim_id','proximal distal'};
    units.summary = {'threshold (uA)','velocity (m/s)'};
    units.spike_counts = 'count, axons x stimulus';
                   
    stimulus.timing = stimulus.t;
    stimulus.pattern = stimulus.p;
    stimulus.current = stimulus.a;
    stimulus = rmfield(stimulus,{'a','t','p'});
    
    notes = { 'axon file:', tools.file('T',[AX.folder AX.name]), ...
              'eidors (v_extracellular) file:', EM.filename, ...
              'output file:', output_file, ...
              'called from:', working_dir };
      
    if ~exist(fileparts(output_file),'dir'), mkdir(fileparts(output_file)); end
    if exist(output_file,'file'), 
      warning('overwrite %s', output_file)
      error TODO_fix_this_tools.file('get',~~~)
    end

    v_list = {'pop','axon_index','diameter','threshold','velocity','g_ratio', ...
              'notes','units','stimulus','results','spike_counts','spike_fraction'};
    
    for v = 1:length(v_list)
      if ~exist(v_list{v},'var'), v_list{v} = ''; end
    end
    
    v_list(cellfun(@isempty, v_list)) = [];
    printf('Saving %s\n',output_file)
    save(output_file,v_list{:});
    clear(v_list{2:end}) % remove variables from workspace
    
  end
end

%% Cleanup

cd(working_dir);

if nargout > 0
  % goal_metric = log10(abs(1-nanmean(goal_metric_vector))); 
  goal_metric = -nanmean(goal_metric_vector); % gradient /descent/
end

return


%% Parse input argumens for nerve_stimulation
function [pop,stimulus] = process_input_args(get_,named,EM,AX)
% Process input arguments, generates nE, nF, fascicle_list, ... 

pop = AX.pop;

nF = sum(strncmp(EM.model.object_name,'Fascicle',4));
nE = numel(EM.model(1).stimulation); 

if any(named('-fasc')), fascicle_list = get_('-fasc');
 if ~any(named('-q')), str_FL = sprintf(', #%d', fascicle_list);
     fprintf('Restrincting simulation to fascicles%s\n', str_FL(2:end))
 end
else fascicle_list = 1:nF; 
end

% If -grid-XY replace axon positions with XY grid
if any(named('-xy')),
 if ~any(named('-q')), disp('using XY grid'), end
 pop = make_axon_xy_grid(EM,pop,AX.nerve); 
end

abs_ax_coord = mean(abs(reshape(abs(cat(1,pop.axon_xy)),[],1))); 
if iscell(EM.info.DomainSize), max_ax_coord = EM.info.DomainSize{1};
else max_ax_coord = EM.info.DomainSize(1);
end

if any(named('-units-um')) || abs_ax_coord > max_ax_coord
  if ~any(named('-units-um'))
    warning('ViNERS:inputAxonUnits', 'Converting output units from �m to mm.')
  end
  pop.fibre_xy = pop.fibre_xy / 1000; % um to mm
  pop.unmyelinated_xy = pop.unmyelinated_xy / 1000; % um to mm
end

if any(named('-downs')), pop = downsample_axons(pop,get_('-downs')); end


%% Set stimulus waveform

has_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if any(named('-stimu')), stimulus = get_('-stimu');
else
  
  stimulus = struct;
  
  pw = 0.1;
  ipg = 0.05;
  
  if any(named('-pw')), pw = get_('-pw'); end
  if any(named('-ipg')), ipg = get_('-ipg'); end
  
  stimulus.t = 30+[0 pw pw+ipg 2*pw+ipg];
  stimulus.p = [1; 0; -1; 0]; 
  stimulus.a = []; % auto
  
  if any(named('-ua')), stimulus.a = get_('-ua'); end
  clear pw ipg idx
end

if ischar(stimulus) % stimulus from file ? 
  if any(stimulus == '~'), stimulus = tools.file(stimulus); end 
  if ~any(named('-q')), fprintf('Loading %s\n', file); end
  if is_ext_(stimulus,'.mat'),      stimulus = load(stimulus); 
  elseif is_ext_(stimulus,'.json'), stimulus = tools.parse_json(stimulus);
  elseif is_ext_(stimulus,'.xml'),  stimulus = tools.parse_xml(stimulus); 
      error TODO_convert_XML_to_struct
  else error('unknown filetype on "%s", expected {.mat, .json, .xml, or <struct>}', stimulus)
  end
end


if any(named('-CL')) % Bionics "CL" equation
  
  bionics_CL_to_uA =  @(CL) 17.5*100.^(CL/255); % convert, units uA
  
  CL = get_('-CL');   
  if ischar(CL), CL = eval(CL); else CL = unique(reshape(CL,[],1)); end
  if numel(CL) == 1, CL = (0:CL:255)'; end % eg -CL 5    
  
  stimulus.a = bionics_CL_to_uA(CL); 
  stimulus.CL = CL;
end


EM_stimPattern = cat(2,EM.model.stimulation.stim_pattern);
isa_monopolar = all(EM_stimPattern(end,:)) || all(EM_stimPattern(end-1,:)); 

cc = 1;
if any(named('-pair')), cc = get_('-pair'); end % vague synonyms
if any(named('-elec')), cc = get_('-elec'); end
if any(named('-chan')), cc = get_('-chan'); end

if any(named('-mono'))
  if ~isa_monopolar, error('%s was generated with bipolar stimuli', EM.filename); end
  assert(size(stimulus.p,2) == numel(cc),'stimulus pattern size must match number of specified channels')
  EM.v_extracellular = EM.v_extracellular(:,cc);
    
elseif isa_monopolar
  if numel(cc) == 1, cc = mod(cc-[1 0],nE)+1; end % default: sequential bipole
  scale = [1; -(ones(numel(cc)-1,1)/(numel(cc)-1))];
  if any(named('-pattern')), scale = reshape(get_('-pattern'),numel(cc),[]); end
  EM.v_extracellular = EM.v_extracellular(:,cc) * scale;   
  if ~any(named('-q')), 
    fprintf('Converting to bipolar stimulus, E%s\n', sprintf('%d',cc))
  end

else 
  if numel(cc) > 1, error('%s was generated with bipolar stimuli', EM.filename); end
  EM.v_extracellular = EM.v_extracellular(:,cc);
end

assert(size(stimulus.p,2) == size(EM.v_extracellular,2), ...
        'stimulus pattern size must match number fields')

%% Put outputs and modified variables in calling workspace 

assignin('caller','nE',nE);
assignin('caller','nF',nF);
assignin('caller','fascicle_list',fascicle_list);
assignin('caller','EM',EM);

return

%% Throw out xy coordinates from axons~/axons.mat and use an XY grid
function pop = make_axon_xy_grid(EM,pop,nerve)

named = evalin('caller','named');
get_ = evalin('caller','get_');

use_F_anatomy = ~any(named('-xy-e')); % XY-em or xy-eidors
if use_F_anatomy, nF = size(nerve.fascicles,3);
else nF = sum(strncmpi(EM.model.object_name,'Fascicle',6));
end

resol = 11; 
use_1d_sample = ~any(named('-xy-g')); % XY-grid or XY-quick
if any(named('-xy-c')), resol = get_('-xy-c'); end  
if any(named('-xy-r')), resol = get_('-xy-r'); end  

%%
if ~use_F_anatomy
  % fac_ =  @(n) EM.(); 
  
  fac_ = @(n) unique( EM.model.elems( ...
                      EM.model.object_id{strcmp(EM.model.object_name, ...
                                        sprintf('Fascicle%d',n))}, : ));
  
  x_ = @(n) EM.model.nodes((n),1); 
  y_ = @(n) EM.model.nodes((n),2); 
  z_ = @(n) EM.model.nodes((n),3); 
end

% source = pop;

for ty = 1:numel(pop)
            
  pop(ty).fibre_diam = median(pop(ty).fibre_diam);
  if pop(ty).myelinated

    if any(named('-ad')), pop(ty).fibre_diam = get_('-ad'); end
    if any(named('-ag')), pop(ty).g_ratio = get_('-ag'); end
    if any(named('-gr')), pop(ty).g_ratio = get_('-gr'); end
      
    pop(ty).g_ratio = median(pop(ty).g_ratio);
    pop(ty).axon_diam = pop(ty).fibre_diam .* pop(ty).g_ratio;
    
    if ~any(named('-q')), 
      fprintf('[%02d | %s] %s, d=%0.4f, g=%0.4f\n', ty, ...
            pop(ty).axon_model,'Simulating grid of axons', ...
            pop(ty).fibre_diam, pop(ty).g_ratio); 
    end        
  elseif ~any(named('-q')), 
      if any(named('-cd')), pop(ty).fibre_diam = get_('-cd'); end
      fprintf('[%02d | %s] %s, d=%0.4f\n', ty, ...
            pop(ty).axon_model,'Simulating grid of axons', ...
            pop(ty).fibre_diam); 
  end
  pop(ty).size_sample = 0; 
  
  if ty > 1      
    pop(ty).axon_xy = pop(1).axon_xy; 
    pop(ty).fascicle = pop(1).fascicle; 
    continue
  end

  %% Get XY positions for each fascicle 
  pop(ty).axon_xy = []; 
  pop(ty).fascicle = []; 
  
  for ff = 1:nF
  
    if use_F_anatomy
      xy_fac = unique(nerve.fascicles(:,:,ff),'rows','stable'); 
    else
      xy_fac = [z_(fac_(ff)) y_(fac_(ff))];
      sel = convhull(xy_fac(:,1),xy_fac(:,2));
      if numel(sel) > 100, sel = sel(round(linspace(1,end,101))); end
      xy_fac = xy_fac(sel([1:end 1]),:);
    end
  
    if numel(resol) > 1 % explicit axon coordinates 

      gx = resol(:,1); 
      gy = resol(:,2); 
      ok = inpolygon(gx(:),gy(:),xy_fac(:,1),xy_fac(:,2));
      
    elseif use_1d_sample % -xy-quick : line sample from min-y to max-y
    
      [~,idx] = min(xy_fac(:,2));
      [~,idx(2)] = max(xy_fac(:,2));
    
      gx = linspace(xy_fac(idx(1),1),xy_fac(idx(2),1),resol+1);
      gy = linspace(xy_fac(idx(1),1),xy_fac(idx(2),1),resol+1);
      gx = conv(gx,[1 1]/2,'valid')'; 
      gy = conv(gy,[1 1]/2,'valid')'; 
      ok = true(size(gx));
    
    else % -xy-grid
    
      xy0 = [median(xy_fac) min(xy_fac) max(xy_fac)];
      xy0(3:4) = max(abs(xy0(1:2)-xy0([3 4; 5 6])),[],1);
      [gx,gy] = meshgrid(xy0(1)-xy0(3)*linspace(-1,1,resol), ...
                         xy0(2)-xy0(4)*linspace(-1,1,resol));
      ok = inpolygon(gx(:),gy(:),xy_fac(:,1),xy_fac(:,2));
    end
  
    pop(ty).axon_xy = [pop(ty).axon_xy; gx(ok) gy(ok)];
    pop(ty).fascicle = [pop(ty).fascicle; 0*gx(ok) + ff];
  end
end
  
n_axons = size(pop(1).fascicle); 
for ty = 1:numel(pop)
  for f = {'fibre_diam','axon_diam','g_ratio'}
    if isempty(pop(ty).(f{1})), continue, end
    pop(ty).(f{1}) = pop(ty).(f{1})(1) * ones(n_axons);
  end
end

return

function [pop] = downsample_axons(pop,fraction)
% [pop,sam] = downsample_axons(pop,sam,get_('-downs'));

nType = numel(pop);
nF = max(cat(1,pop.fascicle));
if numel(fraction) == 1, fraction = repmat(fraction, nType, 1); end

f_sum_ = @(f_id) arrayfun(@(f) sum(f_id == f), 1:nF);
[pop.source_total] = deal(0);
[pop.source_index] = deal(0);

for ty = 1:nType
  
  pop(ty).source_total = f_sum_(pop(ty).fascicle);
  pop(ty).source_index = (1:numel(pop(ty).fascicle))';
  
  if fraction(ty) == 1, continue, end
  if fraction(ty) > 1, fraction(ty) = 1/fraction(ty); end
  
  this = pop(ty);
  
  values = [this.fibre_diam this.axon_diam this.g_ratio this.axon_xy];
  values = zscore(values);
  
  sample = false(size(values(:,1))); 
  
  for ff = 1:nF % get samples within each fascicle 
    
    f_sel = (this.fascicle == ff);
      
    target = ceil(fraction(ty) * sum(f_sel));
    v = values(f_sel,:);
    
    typ = nanmedian(v);
    [~,sel] = min(sum((v-typ).^2,2));  

    extrema = [];
    [~,extrema(:,1)] = min(v); 
    [~,extrema(:,2)] = max(v); 
    extrema = unique(reshape(extrema,[],1)); 

    while numel(sel) < target

      u = arrayfun(@(s) sum((v - v(s,:)).^2,2), sel, 'unif',0);    
      e = arrayfun(@(s) sum((v - v(s,:)).^2,2), extrema, 'unif',0);
      d = mean([u{:} 0.2*[e{:}]],2); d(sel) = nan; 
      [~,ix] = nanmax(d); 
      sel = [sel; ix];  %#ok<AGROW>
    end
    
    f_sel = find(f_sel);
    sample(f_sel(sel)) = true; 
  end
  
  clear u e d v keep sel typ ix st sf target extrema ff f_sel
  
  for y = fieldnames(this)'
    if size(this.(y{1}),1) ~= size(sample,1), continue, end
    this.(y{1}) = this.(y{1})(sample,:);
  end
  
  this.source_index = find(sample);
  
  pop(ty) = this;
end

%% All axon populations with a certain axon model
function [pop,sam] = merge_populations(pop,axon_type)

idx = arrayfun(@(p) p*ones(size(pop(p).size_sample)), 1:numel(pop),'unif',0 );  
[pop.population_id] = deal(idx{:});
sel = strcmp({pop.axon_model},axon_type{1}); 
pop = pop(sel);   

for f = fieldnames(pop)'
  if numel(pop) <= 1, break, end
  if ischar(pop(1).(f{1})), pop(1).(f{1}) = {pop.(f{1})};
  else pop(1).(f{1}) = cat(1,pop.(f{1})); 
  end
end

pop = pop(1); pop(1).axon_model = axon_type{1};
assert(all(pop.myelinated) == any(pop.myelinated),'Myelination is all-or-none')
pop.myelinated = pop.myelinated(1);

if nargout == 1, return, end

% Since each axon gets simulated independently anyway, it's not important
% to merge across myelinated classes, unmyelinated classes for the purposes
% of sampling the axon diameters intelligently (compare merge_populations
% in models.membrane_currents)

sam = struct;
sam.index = pop.size_sample;
sam.count = arrayfun(@(x) sum(sam.index == x), 1:nG)';
the = @(v) arrayfun(@(x) median(pop.(v)(sam.index == x)), (1:nG)');

if pop.myelinated
  sam.fibre_diam = the('fibre_diam');
  sam.axon_diam = the('axon_diam');
  sam.g_ratio = the('g_ratio');
else
  sam.fibre_diam = the('fibre_diam');
end




