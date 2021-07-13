
function get_example_cachefiles ( variables, user_examples, user_indices )      
% This function is part of the standard suite of local result caching tools
% implemented by Calvin Eiber for the MATLAB + NEURON simulation

% models.axon_model generate output files in ./cache/ (or the path defined 
% by tools.cache). This function gathers a selection of examples of these
%  output files which covers the range of simulated behaviours. 

% version 0.3  03-Jul-2020  Calvin Eiber > refactor to use tools.cache

sel_name = {}; 

if nargin == 2, 
  q_val = user_examples(:)';
  q_name = arrayfun(@(q) sprintf('q%0.2f_',q), q_val,'Unif',0);    
  for ii = 1:numel(variables)
    sel_name = [sel_name strcat(q_name, variables{ii})]; %#ok<AGROW>
  end  
else 
  for ii = 1:numel(variables)
    sel_name = [sel_name strcat({'min_','median_','max_'}, variables{ii})]; %#ok<AGROW>
  end  
  q_val = [0 0.5 1];
end

sel_index = ones(size(sel_name)); 

for ii = 1:length(sel_name)    
    if mod(ii,length(q_val)) == 1
        y = evalin('caller',variables{((ii-1)/length(q_val)+1)}); 
    end    
    
    q = q_val(mod(ii-1,length(q_val))+1); 
    [~,sel_index(ii)] = min(abs(y - quantile(y,q)));
end

if nargin > 2, 
    sel_name  = [ sel_name  user_examples ]; 
    sel_index = [ sel_index user_indices  ]; 
end

[ux,~,~] = unique(sel_index); 

example_data = []; 

list = dir(tools.cache('path','*_out.mat')); 
if ~isfield(list,'folder'), % octave compatibility
  [list.folder] = deal(tools.cache('path'));
end

list_id = arrayfun(@(n) str2double(regexp(n.name,'-?\d+','match','once')), list);
f_ = @(f) [f.folder filesep f.name];

if mean(diff(sort(list_id))) > 1
 if evalin('caller','exist(''axon_index'',''var'');')
  axon_id = evalin('caller','axon_index');
  ux = axon_id(ux); 
 else warning('get_cacheFiles:axonIndex', ...
              ['The indices of the cache files look out of sequence. ' ... 
               'Defining "axon_index" in my calling function might help ' ...
               '(%d files, index range %03d - %03d).'], numel(list_id), ...
                    min(list_id), max(list_id))
  if max(ux) <= numel(list_id), 
      axon_id = sort(list_id);
      ux = axon_id(ux); 
  end
 end  
end

for ii = 1:length(ux)
    
    if ~any(list_id == ux(ii)) % requested cache file doesn't exist
        fprintf('[%cMissing file: n%03d_out.mat]%c\n',8,ux(ii),8)
        ux(ii) = -1; 
        continue
    end
    
    this = load(f_(list(list_id == ux(ii))));
    % this.sweep_ID = ux(ii);    % appears to be wrong ???  
    this.result_unit = {'# spikes, no stim','Threshold, uA', ...
                        'Velocity, m/s','Total length, um'};
    if isempty(example_data), example_data = this;
    else                      example_data(end+1) = this; %#ok<AGROW>
    end
end

if any(ux == -1) % some missing files ... 
    sel_name(ux == -1) = []; 
    sel_index(ux == -1) = [];
end

assignin('caller', 'selected_examples', sel_name)
assignin('caller', 'example_index',     sel_index)
assignin('caller', 'example_data',      example_data)