
function status =  model_status(stage, label, varargin)
% tools.model_status( stage, label, ... ) returns the model status 
%   ( 'done', 'incomplete', 'missing', 'not ready' ) for ViNERS for a given 
%   step & output label. The part of the filename in (parens) is the label. 
%   If label is empty, any model output from that stage will do. 
% 
%  The model stages are: 
%    setup:
%       axon_population (axons)
%       nerve_anatomy (sensitivity or potential, as appropriate)
%       membrane_currents (I_m)
%    outputs: 
%       axon_sfap
%       axon_thresholds
%       nerve_stimulation (stim)
%       nerve_recording   (ENG) *
%       ecap_recording    (ecap)
% 
% tools.model_status( stage, label, 'done' ) returns true/false depending 
%  on whether a certain stage has been completed. Simiarly, 
% tools.model_status( stage, label, 'ready' ) returns true/false depending 
%  on whether the prerequisites for a certain stage have been met. 
% 
% Be sure to set tools.file('set','sub-xxx') before running this to make
%  sure you're working in the correct subject, sample, and run! 
% 
% NOTE: for nerve_recording, the syntax is: 
%           model_status(stage, waveType, label, ...)
% see models.nerve_recording('-list') for a list of ENG recording types. 
% 
% v0.1 CDE 24 June 2021

if nargin < 2, label = '';
elseif strcmp(label,'done'), label = ''; varargin = [{'done'} varargin]; 
end

named = @(v) strncmpi(v,varargin,length(v)); 
% get_ = @(v) varargin{find(named(v))+1};

    
if any(named('-q')), printf = @(varargin) []; else printf = @fprintf; end
printf('%s\n', tools.file('T',tools.file('sub~'))); 

model_steps = { {'axons','axon_pop','axon_population'}, ... % model set-up steps
                {'sensitivity','sens','T_e','nerve_anatomy'}, ...
                {'potential','pot','V_e','nerve_anatomy -stim'}, ... 
                {'currents','I_m','mem','membrane_currents'}, ...
                {'sfap','axon_sfap'}, ...
                {'thresholds','axon_t','axon_thresholds'}, ...
                {'stimulation','stim','nerve_s','nerve_stimulation'}, ...
                {'recording','ENG','wave','nerve_r','nerve_recording'}, ...
                {'ecap','ECAP_recording'}}; % , ...
                % {'raster','spikes','random_raster'} };

if 0
  %% labels should be uniquely specifying
  disp('label validation:') %#ok<UNRCH>
  validate = [model_steps{:}]; 
  for ii = 1:numel(validate), 
    if sum(strncmpi(validate,validate{ii},length(validate{ii}))) > 1 
      disp(validate{ii})
    end
  end
end

if nargin == 0, stage = cellfun(@(x) x(1), model_steps(1:end)); end
if ~iscell(stage), stage = {stage}; end % handle multiple stage queries
status = cell(size(stage)); 

%% Get files
list = dir(tools.file('sub~/**/*.mat')); subs = tools.file('sub~/');
list = arrayfun(@(x) [strrep(x.folder,subs,'/') filesep x.name], list, 'unif', 0);
list = strrep(list,filesep,'/'); % unix-style seperators

%% Parse requests 
parse_request(); % reset cache
for step_id = 1 : numel(model_steps)
    
  % request is a vector matching status() and stage()
  request = cellfun(@(a) any(cellfun(@(b) any(strncmpi(a, b, numel(b))), ...
                                   model_steps{step_id})), stage);  
  if ~any(request), continue, end 
  if ismember(step_id, 8) % nerve_recording: special syntax
    if nargin > 2, w_lbl = varargin{1}; else w_lbl = ''; end
    status(request) = parse_request(list,label,step_id,w_lbl);    
  else
    status(request) = parse_request(list,label,step_id);
  end
end

%% Display

for step_id = 1 : numel(model_steps)
  if any(named('-q')), break, end  
  request = cellfun(@(a) any(cellfun(@(b) any(strncmpi(a, b, numel(b))), ...
                                   model_steps{step_id})), stage); 
  if ~any(request), continue, end
  printf('%-20s %s\n', model_steps{step_id}{end}, status{find(request,1)})
end

if any(named('done')),      status = strcmp(status,'done');
elseif any(named('ready')), status = ~strcmp(status,'not ready');
elseif nargout == 0, clear
end

return





function status = parse_request(list,label,index,varargin)

persistent cache
if nargin == 0, get_axons_file(); cache = {}; return, end % reset caches
if numel(cache) >= index && ~isempty(cache{index}), % answer from cache
    status = cache{index}; return, 
end

find_ = @(src,expr) cellfun(@(f) ~isempty(regexp(f,expr,'once')), src);

status = []; 

while 1 % switch index ... break 
 switch index 
  case 1 % axon_population     
    status = simple_check(list, label, '/axons/ax[^/]*.mat');
  case 2 % nerve_anatomy (T_e)
    status = simple_check(list, label, '/eidors/sens[^/]*.mat');
  case 3 % nerve_anatomy -stimulation (V_e)    
    status = simple_check(list, label, '/eidors/stim[^/]*.mat');
  case 4 % membrane_currents
    %%  
    d = get_axons_file(list, label); 
    if isempty(d), status = {'not ready'}; break, end
    if isempty(label), in = '[^/]*'; else in = label; end
     
    axon_models = unique({d.pop.axon_model}); 
    for aa = 1:numel(axon_models) % check axons 
      nG = max(cat(1,d.pop(strcmp({d.pop.axon_model},axon_models{aa})).size_sample));
      n_items = sum(find_(list,['^/axons/' in '/' axon_models{aa} '/[^/]*'])); 
      if n_items == nG+1, continue
      elseif aa == 1 && n_items == 0, status = {'missing'}; 
      else                            status = {'incomplete'};
      end, break
    end
    if isempty(status), status = {'done'}; end
     
   case 5 % axon_sfaps
     %%
     status = simple_check(list, label, '/sfap/sens[^/]*.mat');
     if strcmp(status{1},'missing'), % are prereqs met?
       d = get_axons_file(list, label); 
       if isempty(d), status = {'not ready'}; break,  end
       ok = simple_check(list, '', '/eidors/sens[^/]*.mat');
       if ~strcmp(ok{1},'done'), status = {'not ready'}; end
     end
     
   case 6 % axon_thresholds
     %% Find /thresholds/label folder
     d = get_axons_file(list, label); 
     if isempty(d), status = {'not ready'}; break, end
     if isempty(label), in = '[^/]*'; 
     elseif ~any(label == '('), in = ['[^/]*\(' label '\)[^/]*']; 
     else in = label; 
     end     
     subset = list(find_(list,['^/thre[^/]*/' in]));
     if isempty(subset), 
       ok = simple_check(list, '', '/eidors/stim[^/]*.mat');
       if strcmp(ok{1},'missing'), status = {'not ready'};
       else                        status = {'missing'};
       end,                        break
     end
     %% Check each fascicle 
     nF = size(d.nerve.outline,3); 
     axon_models = unique({d.pop.axon_model}); 
     for aa = 1:numel(axon_models) % check each model
        n_items = sum(find_(subset,axon_models{aa}));
        if n_items ~= nF, 
            status = {sprintf('incomplete (%s %d/%d fascicles)', ...
                               axon_models{aa}, n_items, nF)}; 
            break
        end
     end
     if isempty(status), status = {'done'}; end

   case 7 % nerve_stimulation
       
     d = get_axons_file(list, label); 
     if isempty(d), status = {'not ready'}; break, end

     if isempty(label), in = '[^/]*'; 
     elseif ~any(label == '('), in = ['[^/]*\(' label '\)[^/]*']; 
     else in = label; 
     end     
     subset = list(find_(list,['^/stim[^/]*/' in]));
     
     if isempty(subset), 
       ok = simple_check(list, '', '/eidors/stim[^/]*.mat');
       if strcmp(ok{1},'missing'), status = {'not ready'};
       else                        status = {'missing'};
       end,                        break
     end
     
     %% Check each fascicle
     nF = size(d.nerve.outline,3); 
     axon_models = unique({d.pop.axon_model}); 
     for aa = 1:numel(axon_models) % check each model
        n_items = sum(find_(subset,axon_models{aa}));
        if n_items < nF, 
            status = {sprintf('incomplete (%s %d/%d fascicles)', ...
                               axon_models{aa}, n_items, nF)}; 
            break
        end
     end
     if isempty(status), status = {'done'}; end

   case 8 % nerve_recording
     
     assert(nargin > 3,'not enough input arguments')
     w_type = label; 
     w_lbl = varargin{1}; 
     if isempty(w_type), 
       if isempty(w_lbl), in = '[^/]*'; 
       else in = ['[^\(]*\(' w_lbl '[^/]*']; 
       end
     elseif any(w_type == '('), in = [w_type '[^/]*'];
     elseif isempty(w_lbl), in = [w_type '[^/]*']; 
     else in = [w_type '[^\(]*\(' w_lbl '[^/]*'];        
     end
     
     subset = list(find_(list,['^/waves[^/]*/' in]));
     subset(find_(subset,'^/waves[^/]*/stim')) = []; % models.ecap_record

     if isempty(subset)
       ok = parse_request(list,'',4); % membrane currents?
       if strcmp(ok{1},'missing'), status = {'not ready'}; break, end
       
       ok = simple_check(list, '', '/eidors/stim[^/]*.mat');
       if strcmp(ok{1},'missing'), status = {'not ready'};
       else                        status = {'missing'};
       end,                        break
     end
     
     par = models.nerve_recording('-list');
     sel = cellfun(@(n) contains(subset{1},n), {par.name});
     par = par(sel);
     
     if isempty(par), 
       status = {'incomplete (?? missing from list ??)'}; 
       break
     end % not listed
     
     n_expected = numel(par.spikerate) * numel(par.coherence) * ...
                  numel(par.frequency) * numel(par.exponent) * par.n_reps;
     if numel(subset) == 45 && contains(subset{1},'flat'), 
       status = {'done'}; % default produces flat/(45 files)
     elseif numel(subset) < n_expected, 
       status = {sprintf('incomplete (%d/%d waves)', ...
                          numel(subset), n_expected)}; 
     else status = {'done'};
     end
     
     
   case 9 % ECAP_recording
              
     if isempty(label), in = '[^/]*'; 
     elseif ~any(label == '('), in = ['[^/]*\(' label '\)[^/]*']; 
     else in = label; 
     end     
  
     subset = list(find_(list,['^/waves/stim[^\(]*' in]));

     if isempty(subset)
       ok = parse_request(list,'',7); % nerve stim?
       if strcmp(ok{1},'missing'), status = {'not ready'}; break, end
       
       ok = parse_request(list,'',4); % membrane currents?
       if strcmp(ok{1},'missing'), status = {'not ready'}; break, end

       ok = simple_check(list, '', '/eidors/sens[^/]*.mat');
       if strcmp(ok{1},'missing'), status = {'not ready'};
       else                        status = {'missing'};
       end,                        break
     end

     % Check how many waves there are
     d = load([tools.file('sub~') subset{1}]);
     
     if numel(subset) < numel(d.stimulus.current)
          status = {sprintf('incomplete (%d/%d waves)', ...
                          numel(subset), numel(d.stimulus.current))}; 
     else status = {'done'};
     end
     
   case 10 % random_raster   
     warning('Not implemented yet')
     
   otherwise error('case %d not implemented', step_id)
 end
 break
end

cache{index} = status; 









function [status,file] = simple_check(list, label, expr) 

 find_ = @(src,e) cellfun(@(f) ~isempty(regexp(f,e,'once')), src);

 file = list(find_(list,['^' expr]));
 if isempty(file), status = {'missing'}; return, end
 if isempty(label), status = {'done'};  return, end
 file = file(find_(file,['\(' label '\)']));

 if isempty(file), status = {'missing'};
 else             status = {'done'}; 
 end

function data = get_axons_file(list, label)
     
persistent d
if nargin == 0, d = []; return, end
if isempty(d)

    [~,d_name] = simple_check(list, label, '/axons/ax[^/]*.mat');
    if numel(d_name) > 1
        tfs = tools.file('sub~/');
        fileinfo = cellfun(@(x) dir([tfs x]), d_name); 
        [~,sel] = max([fileinfo.datenum]);
        d_name = d_name(sel); 
    end
    
    if isempty(d_name), d = []; 
    else d = load([tools.file('sub~/') d_name{1}]);
    end
end
data = d;

