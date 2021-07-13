
function result = axon_model(index, field, varargin)
% Run multithreaded NEURON simulation
% the expectation is that the model is called in a parfor loop like so: 
% 
% parfor ii = 1:nSimulations 
%   result(ii,:) = axon_model(ii, VE_field or [],arg_in{:});
% end
% 
% Possible command-line arguments: 
% -get-parameters : return list of biophyics parameters (matches inputs
%        needed by makeFromTemplate for +models/Gaines2016_v3.hoc.template)
%        these parametes can be set using 'Name',value semantics. 
% -sensitivity  [values] : accept varlist for sensitivity analysis of
%                                  biophysical parameters (Sobol')
%
% -sphase      : run spatial phase sensitivity test
% -stimulus    : run a specified current level (instead of threshold)
% -constantlen : set length of model axon to constant value (recommended)
% -gRatio      : set axon g-ratio (recommended, irrelevant to -sundt)
% 
% -motor; MRG  : set model type to (modified) MRG motor axon model
% -sensory     : set model type to Gaines sensory axon model
% -cfiber      : set model type to Sundt nociceptor model
% -XY          : coordinates of axon within Ve field
% 
% -sweep-range : sweep range for find_threshold_external / _internal
%                (1x3 vector)
% 
% This depends on C:\neuron\bin\nrniv.exe existing and that nrniv.exe can 
% be executed from from the system command-line. Executing the neuron model
% can be slow; however, this code will hang in near line 96 if there are
% any problems executing the neuron model. If this code takes excessively
% long, try running with -debug to view the neuron output (the system call
% which invokes neuron is run with -echo on)
% 
% model version 0.4 Calvin Eiber 03-Jun-2020

if nargin == 0, index = 0; field = []; varargin{1} = '-spike';  end
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

contains = @(a,b) ~isempty(strfind(a,b)); %#ok<STREMP>

if ~contains(ctfroot, 'MATLAB')
  % [~,pidlist] = system('pgrep octave'); 
  % pidlist = str2double(regexp(pidlist,'\d+','match'));
  t.ID = getpid; % sum(pidlist < getpid);
  more off % disable octave paging
  warning off Octave:regexp-lookbehind-limit  
else t = getCurrentTask(); % Thread index 
  if isempty(t), t.ID = 0; end
end

NOT_MYEL = any(named('-unmy')) || any(named('Sundt')); 
FLAG_DEBUG = any(named('-debug')); 

% Get list of biophysics parameters and parse input args
vars = get_paramaters(NOT_MYEL, t.ID);

if any(named('-get-par')), result = vars'; return, end

%% Build files and template

cache_path = tools.cache('path'); 
cache_file = sprintf('%s/n%05d_out.mat',cache_path,index);
hoc_file   = sprintf('%s/nrn-%02d.hoc',cache_path,t.ID); 

% the Sundt model does /not/ include the KCNQ channels except at
% the T-junction in the DRG. 

if NOT_MYEL, template_file = {'./+models/Sundt2015_v1.hoc.template'};
elseif any(named('-sens')) || any(named('gaines'))
     template_file = {'./+models/Gaines2016_v3.hoc.template', 'neuron_SENSORY'};
else template_file = {'./+models/Gaines2016_v3.hoc.template', 'neuron_MOTOR'};
  if ~any(named('mrg')), 
    warning('ViNERS:axonModel:default','using default axon model (modified from MRG 2002)')
  end
end

cache_path = regexprep(cache_path,'[/\\]','/'); % makeFomTemplate garbles the other slash, and windows doesn't care
tools.makeFromTemplate(template_file{:},vars{:},'CACHEPATH',cache_path, ... 
                                                '-output',hoc_file); 
%%
node = make_Vfield_datfiles(t.ID,field,vars,1);

if any(named('-echo')) || FLAG_DEBUG, do_echo = '-echo'; 
else do_echo = false; 
end

if any(named('-sweep-range')), do_range = [1;1] * get_('-sweep-range');
elseif NOT_MYEL, do_range = [0.01 0 0.1; 1000 0 2000];
else             do_range = [0.05 0 0.6;  100 0  400]; 
end

if any(named('-VCLAMP'))
  hoc_cmd = sprintf('run_VoltageClamp()');
  hoc_cmd = set_named_mode(get_,named,hoc_cmd); % Use a specific output mode? 

elseif isempty(field), % Internal step current stimulus
  hoc_cmd = sprintf('find_threshold_internal(%0.4f,%0.4f,%0.4f)',do_range(1,:));
elseif any(named('-fixed-s'))
  hoc_cmd = sprintf('run_stim_sweep()');
  hoc_cmd = set_named_mode(get_,named,hoc_cmd); % Use a specific output mode? 
else
  hoc_cmd = sprintf('find_threshold_external(%0.4f,%0.4f,%0.4f)',do_range(2,:));
end

% ======================================================

nrn_out = invoke_neuron(hoc_file,hoc_cmd,do_echo);

% ======================================================

if any(named('-no-analysis')), result = nrn_out; return, end

% for THRESHOLD_HIGH, check 5x inital range as well
if contains(nrn_out,'abort THRESHOLD_HIGH') 
  if isempty(field)
       hoc_cmd = sprintf('find_threshold_internal(%0.4f,%0.4f,%0.4f)',do_range(1,:)*5);
  else hoc_cmd = sprintf('find_threshold_external(%0.4f,%0.4f,%0.4f)',do_range(2,:)*5);
  end
  nrn_out = invoke_neuron(hoc_file,hoc_cmd,do_echo);
end

[result,hoc_cmd,do_post_analysis] = parse_NRN_threshold(hoc_cmd, nrn_out);

if do_post_analysis % model completed successfully 

  u = invoke_neuron(hoc_file,hoc_cmd,do_echo); %#ok<NASGU>
  [waves, cv] = compute_conduction_velocity(t.ID,node); 
  if ischar(result), result = [compute_apparent_threshold(waves,vars) cv];
  else result(2) = cv;
  end
  
  if any(named('-postprocess')) % custom postprocesser
    custom_cache_file = get_('-postprocess');
    [waves,result] = custom_cache_file(cache_file,waves,vars,result); 
  end
  
  if ~isempty(waves) % generate "default" cache file
    result(3) = range(waves.z);
    if any(named('-term')), waves.prop_segment = numel(node.NODE)-1;
    end
    [waves,result] = build_cache_file(cache_file,waves,vars,result);
  end  
else waves.spiketimes = []; % for struct output  
end
  
if ~any(named('-r-v')) && ~isstruct(result)
  
  r_vals = result; 
  result = struct; 
  result.threshold = r_vals(1);
  result.velocity  = r_vals(2);
  result.spiketime = waves.spiketimes; 
end

return

%% Core function to invoke nrniv.exe 
function output = invoke_neuron(hoc_file,command,echo)
  if nargin == 1, command = 'quit()'; end
  if isempty(command), output = ''; return, end
  
  % C:\neuron\bin\nrniv.exe
  persistent nrn_exe
  if isempty(nrn_exe), nrn_exe = tools.configuration('neuron'); end
  command = sprintf('%s -nobanner -nopython "%s" -c "%s"', nrn_exe, ...
                                          hoc_file, command);
  if nargin > 2 && any(echo), system([command '&']); end
  [~,output] = system(command);
return

%% Parse inputs to models.axon_model
function list = get_paramaters(NOT_MYEL,threadID)

named = evalin('caller','named');
get_ = evalin('caller','get_');
verbose = any(named('-debug')); 

if NOT_MYEL

    list = {'thread','fibreDiam','nodeLength','nNodes','simDuration', ...
            'p_rhoA','p_gKCNQ','p_gNaV','p_gKV','p_NaVshift'};
else
    list = {'thread','fibreDiam','numbernodes','simDuration', ...
            'p_nodelen','p_mysalen','p_flutlen','p_stinlen', ...
            'p_nodedia','p_mysadia','p_flutdia','p_stindia',...
                        'p_mysapax','p_flutpax','p_stinpax','p_nlamella'};
end            
[list{2,:}] = deal(1);

%% Parse input args into structure
list{2,1} = threadID; % thread

if NOT_MYEL
  list{2,2} = 0.8;  % fibreDiam, um
  list{2,3} = 20;   % nodeLength, um
  list{2,4} = 300;  % n_nodes
  list{2,5} = 100;  % simDuration, ms
  if any(named('-old-resol'))
    list{2,3} = 10;  
    list{2,4} = 600; 
  end
else
  list{2,2} = 4;    % fibreDiam, um
  list{2,3} = 31;   % numbernodes
  list{2,4} = 100;  % simDuration, ms
end

for ii = 1:size(list,2) % 'name', value semantics
  if any(named(list{1,ii})), list{2,ii} = get_(list{1,ii}); 
    if verbose, fprintf('%s = %0.6f\n', list{:,ii}), end    
  end
end

if any(named('-sensitivity'))
  param = get_('-sensitivity');
  param(end+1 : size(list,1)) = 1; 
  if NOT_MYEL, list(2,6:9)  = num2cell(param(:));
  else                list(2,5:16) = num2cell(param(:)); 
  end
end

if any(named('-constantlen'))
  param = get_('-constantlen');
  if NOT_MYEL
    list{2,4} = ceil(1000 * param / list{2,3}) + 1;
    if verbose, fprintf('%s = %d\n', list{:,4}); end
  else 
    node = make_Vfield_datfiles(threadID,[],list,0);
    list{2,3} = round(list{2,3} * param / range(node.NODE));
    if verbose, fprintf('%s = %d\n', list{:,3}); end
  end
end

if any(named('g_ratio')) && ~NOT_MYEL
    list = set_gRatio_scale(list, get_('g_ratio'));
elseif any(named('-gra'))
    list = set_gRatio_scale(list, get_('-gra'));
end

if any(named('-sphase')) % spatial phase sensitivity test
  param = get_('-sphase');
  list{1,end+1} = 'spatial_phase';
  list{2,end} = param; 
end

if any(named('-XY')) % coordinates of axon within field
  param = get_('-XY');
  list{1,end+1} = 'spatial_location';
  list{2,end} = param;
  if size(param,2) == 3 || isstruct(param), list = convert_3D_axon(param,list); end
end

if any(named('-stim')) % approach for handling complex stimulation
  param = get_('-stim');
else % construct a biphasic 1-channel stimulus
  
  delay = 30; 
  p_width = 0.4;
  
  if any(named('-VCLAMP'))    
    delay = 50; 
    p_width = 20;
  end
  
  if any(named('-pw')), p_width = in_('-pw'); end
  if any(named('-delay')), delay = in_('-delay'); end


  param = struct;
  param.t = delay + [0 p_width p_width+0.02 2*p_width+0.02];
  param.p = [1; 0; -1; 0]; 
  param.a = unique(reshape([1;10]*round(logspace(1,2,13),1),[],1));
  
  if any(named('-ua')), param.a = get_('-ua');
  elseif any(named('-VCLAMP'))
    if any(named('-inactivation')), param.a = -120:5:0;
    else                            param.a = -80:5:60;
    end  
  else    
    param.a = unique(reshape([1;10]*round(logspace(1,2,13),1),[],1));
  end

end
list{1,end+1} = 'stimulus_pattern';
list{2,end} = param; 
function vars = set_gRatio_scale(vars, gratio)

% nodeD  = (0.00630378*(fiberDiam*fiberDiam) + 0.207054*(fiberDiam) + 0.5339) * $p_nodedia // quadratic re-scaling per Lubba 2019
% axonD  = (0.01876226*(fiberDiam*fiberDiam) + 0.478749*(fiberDiam) + 0.1204) * $p_stindia
% paraD1 = (0.00630378*(fiberDiam*fiberDiam) + 0.207054*(fiberDiam) + 0.5339) * $p_mysadia
% paraD2 = (0.01876226*(fiberDiam*fiberDiam) + 0.478749*(fiberDiam) + 0.1204) * $p_flutdia

diam  = vars{2,2}; % fibreDiam; 

axonDiamFcn = @(d) (0.01876226*(d*d) + 0.478749*(d) + 0.1204); 
p_stindia = gratio * diam / axonDiamFcn(diam);

vars{2,9}  = p_stindia;  % p_nodedia
vars{2,10} = p_stindia;  % p_mysadia
vars{2,11} = p_stindia;  % p_flutdia
vars{2,12} = p_stindia;  % p_stindia 

return
function cmd = set_named_mode(get_,named,cmd)

if any(named('-VClamp-full'))
    cmd = ['output_mode=6" -c "' cmd];    
elseif any(named('-spike')) || any(named('-save-S')) 
    cmd = ['output_mode=1" -c "' cmd];
elseif any(named('-volt')) || any(named('-save-V')) 
    cmd = ['output_mode=2" -c "' cmd];
elseif any(named('-current')) || any(named('-save-I')) 
    cmd = ['output_mode=3" -c "' cmd];
end

if any(named('-t-start'))  
    cmd = sprintf('start_output=%0.4f" -c "%s', get_('-t-start'), cmd);
end
if any(named('-spk-volt'))  
    cmd = sprintf('AP_threshold=%0.4f" -c "%s', get_('-spike-volt'), cmd);
end

if any(named('-VCLAMP')) % pass in extra commands 
  
% VC_nodeID = 1 // which node is VClamp attached to? (settable)
% VC_parID  = 1 // which VClamp.amp[0,1,2] is updated? 
  
  if any(named('-VClamp-node'))
    cmd = sprintf('VC_nodeID=%d" -c "%s', get_('-VClamp-node'), cmd);    
  end
  if any(named('-inactivation'))    
    cmd = ['VC_parID=0" -c "' cmd];
  end  
end

cmd = [cmd '" -c "quit()']; % make sure neuron exits


function param = convert_3D_axon(xyz,param)

NOT_MYEL = evalin('caller','NOT_MYEL'); 
par_ = @(s) strcmpi(param(1,:),s);

if NOT_MYEL

  s_len  = param{2,par_('nodeLength')};
  n_node = param{2,par_('nNodes')};
  s_type = 1; 
  
else
  threadID = evalin('caller','threadID');
  ultrastructure = make_Vfield_datfiles(threadID,[],param,'-info');

  s_len = ultrastructure.lengths;  
  s_type = ultrastructure.types;
  n_node = ultrastructure.n_node;
end

s_type = [repmat(s_type,1,(n_node-1)) 1]; % type IDs of each segement
z = cumsum(s_len(s_type))' / 1000;

if isstruct(xyz) && isfield(xyz,'xyz') && ~iscell(xyz.xyz), xyz = xyz.xyz;    
elseif isstruct(xyz)
  
  dz = mean(diff(z));
  axon = xyz; 
  
  eqs_xyz = []; 
  b_index = []; 
  
  for ii = 1:numel(axon.xyz)
    
    if size(axon.xyz{ii},1) > 1
      u = cumsum([0; sqrt(sum((axon.xyz{ii}(1:end-1,:) - axon.xyz{ii}(2:end,:)).^2,2))]); 
      if any(diff(u) == 0)  
        cull = [false; diff(u) == 0];
        axon.xyz{ii}(cull,:) = []; 
        u(cull) = []; 
      end
      axon.xyz{ii} = interp1(u,axon.xyz{ii},0:dz:max(u));        
    end
    
    bi = size(eqs_xyz,1) + [1 size(axon.xyz{ii},1)];
    eqs_xyz = [eqs_xyz; axon.xyz{ii}(:,1:3)]; %#ok<AGROW>
    b_index = [b_index; bi]; %#ok<AGROW>
  end
  
  idx = sub2ind(size(b_index),axon.list(:,1),axon.list(:,2)+1); 
  axon.list(:,1) = b_index(idx); 
  
  idx = sub2ind(size(b_index),axon.list(:,3),axon.list(:,4)+1); 
  axon.list(:,3) = b_index(idx); 
 
  models.branching_axons('-setup','-axon',axon)
  
  if NOT_MYEL
    param{2,par_('spatial_location')} = eqs_xyz;
    param{2,par_('nNodes')} = size(eqs_xyz,1);
  else
    error TODO_3d_branching_structure_myelinatedAxon_output
  end

  param(:,end+1) = {'branching_axons',axon}; % #IFDEF arg
  return
end

if size(xyz,1) == 1, xyz = permute(xyz,[3 2 1]); end

u = cumsum([0; sqrt(sum((xyz(1:end-1,:) - xyz(2:end,:)).^2,2))]); 

if any(diff(u) == 0)  
  cull = [false; diff(u) == 0];
  xyz(cull,:) = []; 
  u(cull) = []; 
end

dz = s_len(s_type)/2e3; % adjust sample points to middle of segment (currently end-of-segment)
if any(par_('spatial_phase'))
 if ischar(param{2,par_('spatial_phase')}), rsp = 2*rand-1;
 else rsp = param{2,par_('spatial_phase')};
 end, dz = dz + diff(z(find(s_type==1,2)))*rsp;
end

eqs_xyz = interp1(u-mean(u),xyz,z-mean(z)-dz','pchip'); 
eqs_xyz(isnan(eqs_xyz(:,1)),:) = [];

param{2,par_('spatial_location')} = eqs_xyz;
if NOT_MYEL, param{2,par_('nNodes')} = size(eqs_xyz,1); end

%%
return
%% Visualisation of results 

clf, hold on
plot3(eqs_xyz(:,1),eqs_xyz(:,2),eqs_xyz(:,3),'-','Color',[0 0 0 0.5],'Clipping','off')
scatter3(eqs_xyz(:,1),eqs_xyz(:,2),eqs_xyz(:,3),[],s_type,'s','filled','Clipping','off')
axis equal, grid on



%% interpolate Ve to get the extracellular (stimulation) potentials
function coords = make_Vfield_datfiles(tid,field,param,iStim)

NOT_MYEL = evalin('caller','NOT_MYEL'); 
par_ = @(s) strcmpi(param(1,:),s);

if nargin < 3, iStim = 1; end % global scaling factor now that iterating iStim up/dn 

if NOT_MYEL, node = {'NODE'};
    
    s_len(1) = param{2,3};
    n_node   = param{2,4};
    s_type   = 1; 
    
else % Myelinated MRG-like double-cable model
    
    node = {'NODE','MYSA','FLUT','STIN'};

    f_diam = param{2,2}; 
    n_node = param{2,3};
    
    s_len(1) = 1.0 * param{2,5};       % node width, um
    s_len(2) = 3.0 * param{2,6};       % mysa width, um
    s_len(3) = (2.5811*(f_diam)+19.59) *param{2,7}; % flut width, um
    s_len(4) = ((221.1322 * f_diam .^ 3.133103) ./ ...
                (f_diam .^ 3.133103 + 593.404881)); % stin width, um
    
    s_type = [1 2 3 4 4 4 4 4 4 3 2];
    if ischar(iStim) && strcmpi(iStim,'-info')
      coords.f_diam = f_diam;
      coords.lengths = s_len; 
      coords.types = s_type;
      coords.n_node = n_node;
      return
    end
end
    
widths = s_len(s_type); 
s_type = [repmat(s_type,1,(n_node-1)) 1]; % type IDs of each segement

z = cumsum(s_len(s_type))' / 1000;
z = z - s_len(s_type)'/2e3; % adjust to centre of each 
z = z - mean(z); % and centre

if any(par_('spatial_phase'))
 if ischar(param{2,par_('spatial_phase')}), rsp = 2*rand-1;
 else rsp = param{2,par_('spatial_phase')};
 end, z = z + diff(z(find(s_type==1,2)))*rsp;
end

if isempty(field), vStim = 0*z; 
else
  
  if any(par_('spatial_location'))
    dz = par_('spatial_location'); 
    
    if size(param{2,dz},2) == 3
      x = param{2,dz}(:,1) ; % / 1000;
      y = param{2,dz}(:,2) ; % / 1000;
      z = param{2,dz}(:,3) ; % / 1000;
    else
      x = 0*z + param{2,dz}(1) ; % / 1000;
      y = 0*z + param{2,dz}(2) ; % / 1000;
    end
  else x = 0*z; y = 0*z; % default: 0,0
  end
  
  if iscell(field), 
    
    vStim = cellfun(@(f) f(x,y,z) * iStim, field, 'unif',0);
    vStim = [vStim{:}];
    
  else vStim = field(x,y,z) * iStim; % extracellular potential
  end
end

if any(isnan(vStim(:)))
    error FIX_NANs % probably a units mismatch
end

folder = tools.cache('path'); 
coords = struct; 

for nn = 1:numel(node), coords.(node{nn}) = z(s_type == nn); end

% if isempty(tid), return, end %% just get coords without saving

if any(par_('stimulus_pattern'))
  
  stim = param{2,par_('stimulus_pattern')};
  time = stim.t;
  save(sprintf('%s/stimulus-time_%d.dat',folder,tid),'time','-ascii')
  
  levels = stim.a;
  if size(levels,2) > 1, s = size(levels); 
       save(sprintf('%s/stimulus-levels_%d.dat',folder,tid),'s','levels','-ascii')  
  else save(sprintf('%s/stimulus-levels_%d.dat',folder,tid),'levels','-ascii')
  end
  
  pattern = stim.p; s = size(pattern);
  save(sprintf('%s/stimulus-pattern_%d.dat',folder,tid),'s','pattern','-ascii')
  
  for nn = 1:numel(node)
    V = vStim(s_type == nn,:); s = size(V);
    save(sprintf('%s/type%d_Ve_%d.dat',folder,nn,tid),'s','V','-ascii')
  end
  
else
  for nn = 1:numel(node)
    return % Code for baseline model pre-Nov-2020 below
    V = vStim(s_type == nn); %#ok<UNRCH>
    save(sprintf('%s/v%s%d.dat',cache,lower(node{nn}),tid),'V','-ascii')
  end
end

%% Analysis Functions for NRN output
function [result,hoc_cmd,followup] = parse_NRN_threshold(hoc_cmd, nrn_out)

field     = evalin('caller','field');
thread_id = evalin('caller','t.ID');
index     = evalin('caller','index');
named     = evalin('caller','named');
get_      = evalin('caller','get_');
vars      = evalin('caller','vars');

followup = false; 

if contains(hoc_cmd,'find_threshold_internal'), noun = 'Spike at';
elseif contains(hoc_cmd,'find_threshold_external'), noun = 'Threshold';
elseif contains(hoc_cmd,'run_stim_sweep') || contains(hoc_cmd,'stimulus')
  
    par_ = @(s) strcmpi(vars(1,:),s);

    stim = vars{2,par_('stimulus_pattern')};
    nStim_init = numel(stim.a);
    nStim_out = numel(strfind(nrn_out,'stimulus[0] ='));
    stim_label = 'Stimulus Sweep ('; 
    
    if all(size(stim.a) > 1)
      nStim_init = size(stim.a,2);
      stim_label = sprintf('Multi Stimulus Sweep (%d pulses / stimulus, ', size(stim.a,1));
    end
    
    if nStim_init ~= nStim_out
      fprintf('[%c[%03d] warning: %s%d requested, %d delivered)]%c\n', ...
              char(8), index, stim_label, nStim_init, nStim_out, char(8))
      log_error(thread_id,index,'stimulus_mismatch',vars,nrn_out);    
    else
      fprintf('[%03d] %s%d stimuli, %0.1f - %0.1f uA)\n', ...
              index, stim_label, nStim_init, min(stim.a(:)), max(stim.a(:)))
    end
    
    followup = true; 
    hoc_cmd = ''; 
    result = nrn_out;
    return
elseif contains(hoc_cmd,'find_threshold'), noun = 'Threshold [OLD]';
else warning('ViNERS:cantParseNRN','I won''t know how to parse the output of -c "%s"', ...
                  hoc_cmd)
    followup = false; 
    hoc_cmd = ''; 
    result = nrn_out;
    return
end


if contains(nrn_out,'threshold=')
    
    followup = true; 
    result(1) = str2double(regexp(nrn_out,'(?<=threshold=)\d+(\.\d+)?','match','once'));        
    fprintf('[%03d] %s: %0.2f uA\n', index, noun, result(1))
    
    if isempty(field)
         hoc_cmd = sprintf('run_stim(0.0, %0.6f)', result(1) * 1.05);
    else hoc_cmd = sprintf('run_stim(%0.6f, 0.0)', result(1) + 1);
    end
    
elseif contains(nrn_out,'abort THRESHOLD_LOW')
    
    log_error(thread_id,index,'low_threshold',vars,nrn_out);    
    result(1) = str2double(regexp(nrn_out,'(?<=LOW at )\d+(\.\d+)?','match','once'));    
    fprintf('[%c[%03d] warning THRESHOLD_LOW, using %0.8f uA]%c\n', char(8), index, result(1), char(8))
    
    if isempty(field),
         hoc_cmd = sprintf('run_stim(0.0, %0.6f)', result(1) * 2);
    else hoc_cmd = sprintf('run_stim(%0.6f, 0.0)', result(1) * 2); 
    end    
    result(1) = 0; 
    followup = true; 
    
elseif contains(nrn_out,'abort THRESHOLD_HIGH')
    result(1:2) = [Inf NaN]; 
    log_error(thread_id,index,'high_threshold',vars,nrn_out);
    fprintf('[%c[%03d] warning THRESHOLD_HIGH ]%c\n', char(8), index, char(8))
else
    log_error(thread_id,index,'bad_model',vars,nrn_out);
    fprintf('[%c[%03d] model failure]%c\n', char(8), index, char(8))
    result(1:3) = NaN;
end

hoc_cmd = set_named_mode(get_,named,hoc_cmd); % Use a specific output mode? 

return
function [d,cv] = compute_conduction_velocity(tid,coords)

filename = sprintf('%s/NEURON_vm%d.dat',tools.cache('path'),tid);
d = tools.read_NRN_datfile(filename); 
is_node = strncmpi(d.secname,'node[',5); 

z_mm = coords.NODE; 

if isfield(coords,'MYSA') % MRG-like model 
     d.z = [coords.NODE; coords.MYSA; coords.FLUT; coords.STIN];
else d.z =  coords.NODE; 
     is_node(:) = true; 
end

t0 = nan(size(d.secname(is_node))); % time of first spike at each node

if ~isfield(d,'spk') && isempty(d.t) % no spikes but spiketimes requested
  vars = evalin('caller','vars');
  par_ = @(s) strcmpi(vars(1,:),s);
  n_stim = size( vars{2,par_('stimulus_pattern')}.a, 1 );
  d.spk = cell(1,n_stim); 
end

if isfield(d,'spk') % spiketime output
  for pp = 1:numel(d.spk) % possibly multiple stimulus levels 
    
    t0 = nan(size(d.secname(is_node))); % time of first spike at each node
    if isempty(d.spk{pp}), continue, end
    
    for seg = 1:numel(coords.NODE)
      
      ts = find(d.spk{pp}(:,2) == seg, 1); 
      if isempty(ts), continue, end
      t0(seg) = d.spk{pp}(ts,1); 
    end
    
    ok = ~isnan(t0);
    ok = ok & diff([-1;t0]) > 0; % orthodromic only

    if mean(ok) > 0.4 && sum(ok) > 3, break, end 
    % grab first stimulus with a decent spike
  end
else vPk = -30; % min(0,median(quantile(d.vm,0.99))); % Vm output
  for seg = 1:numel(coords.NODE) % get time of first spike 
      ts = find(d.vm(:,seg) >=vPk, 1); 
      if isempty(ts), continue, end
      t0(seg) = d.t(ts); 
  end
end

[ts,s0] = min(t0); 
x = abs(z_mm(:) - z_mm(s0)); 
ok = ~isnan(t0);
ok = ok & diff([-1;t0]) > 0; % orthodromic only

if sum(ok) < 2, 
  
  fprintf('[%c[%03d] check conduction velocity!]%c\n', char(8), evalin('caller','index'),char(8))
  
  ok = ~isnan(t0);  
  if sign(nanmean(diff(t0(ok)))) < 1 
    ok = ok & diff([inf;t0]) <= 0; % antidromic only
    fprintf('[%c[%03d] using antidromic spike]%c\n', char(8), evalin('caller','index'),char(8))
  else
    ok = ok & diff([-1;t0]) >= 0; % orthodromic only
  end
  if sum(ok) < 2,   ok = ~isnan(t0);  end  
end

cv = nan; 
if sum(ok) > 1

  if isempty(which('fit'))
    [cv_fit] = ((t0(ok)-ts)*[1 0] + [0 1]) \ x(ok);
    cv = abs(cv_fit(1));
    % this fits t(...) ~= (x(ok) - cv(2)) / cv(1) 
    % TODO confirm this in MATLAB not OCTAVE
  else % use fit('poly1')
    cv_fit = fit(t0(ok)-ts,x(ok),'poly1'); % mm per ms = m per s
    cv = abs(cv_fit.p1); 
  end
  if cv < 0, cv = nan; end
end

d.init_segment = s0; % segment ID of earliest spike location

[~,s1] = max(abs(z_mm - z_mm(s0)));
z1 = mean(z_mm([s0 s1])); 
[~,s1] = min(abs(z_mm - z1)); 
d.prop_segment = s1; % segment ID of "the other direction" 

function threshold = compute_apparent_threshold(d,vars)

has_spk = ~cellfun(@isempty,d.spk);
s0 = find(has_spk,1);

if ~any(s0), threshold = inf; return
else
  
  par_ = @(s) strcmpi(vars(1,:),s);
  
  stimulus = vars{2,par_('stimulus_pattern')};
  delay = min(stimulus.t);
  
  rsl = (s0:numel(d.spk)); % using "rsl" as my scratch variable
  rsl(~has_spk(rsl)) = []; % exclude cases with 0 spikes from list
  rsl = cat(1,d.spk{rsl}) - delay; % gather spiketimes and locations
  
  t0 = min(rsl(rsl(:,1) > 0,1)); 
  t1 = min(d.spk{s0}(d.spk{s0}(:,1) > delay,1)-delay);
  
  rsl = t0 ./ t1; % ratio of spike latencie

  if s0 == 1, threshold = [rsl 1-rsl] * [0; stimulus.a(s0)];
  else        threshold = [rsl 1-rsl] * reshape(stimulus.a(s0-[1 0]),2,1);
  end
end

function [f,result] = build_cache_file(cache_file,waves,vars,result)

NOT_MYEL = evalin('caller','NOT_MYEL'); 
idx = [waves.init_segment waves.prop_segment]; 
par_ = @(s) strcmpi(vars(1,:),s);

f = struct; 

if any(par_('stimulus_pattern'))
  stimuli = vars{2,par_('stimulus_pattern')};
  delay = min(stimuli.t(stimuli.t>0));
else
  stimuli = 'nrn-model-default';
  delay = 0;
end

if isempty(waves.t)
  waves.t = [0 vars{2,par_('simDuration')}]; 
end

if isfield(waves,'im')
  if ~NOT_MYEL
      idx = [idx(1)  idx(2)*[1 2 2 2 2 6 6 6 6 6 6] - ...
                            [0 2 1 2 1 6 5 4 3 2 1]];  
      idx(3:4) = idx(3:4) + find(strncmpi(waves.secname,'MYSA',4),1);
      idx(5:6) = idx(5:6) + find(strncmpi(waves.secname,'FLUT',4),1);
      idx(7:end) = idx(7:end) + find(strncmpi(waves.secname,'STIN',4),1);
  end
end

%% Reformat output structure with more verbose names

f.result = result;
f.parameters  = vars';

f.stimuli = stimuli; 
if isstruct(f.stimuli), f.stimuli.t = stimuli.t - delay; end

f.time = waves.t - delay; 
f.time_unit = 'ms'; 
f.length = waves.z;
f.length_unit = 'mm'; 

f.segment = waves.secname(idx); 
f.index = idx;

f.dt_dx = (mean(diff(f.length(1:idx(2)))) / f.result(2)); % delta-ms per sample

if isfield(waves,'im') % membrane current profiles 
  f.I_membrane = waves.im(:,idx); 
  f.I_unit = 'nA'; % ... (mA/cm�)(�m�) = 10^-11 A, /100 = nA
  f.shift = '@(f,t,n) interp1(f.time,f.I_membrane(:,2:end),t - n*f.dt_dx)';
end

if isfield(waves,'vm') && ~isempty(waves.vm)
  f.V_example = waves.vm(:,idx); 
  f.V_unit = 'mV';
end
  
if isfield(waves,'spk'), d = waves; % alias
  
  f.spiketimes = []; 
  
  d.spk(cellfun(@isempty,d.spk)) = {zeros(0,2)};
  
  f.spikes.time = cellfun(@(s) s(:,1)-delay, d.spk, 'unif', 0);
  f.spikes.node = cellfun(@(s) s(:,2), d.spk, 'unif', 0);
  f.spikes.init = cell(size(d.spk)); 
  
  f.spiketimes = cellfun(@(s) s(s(:,2) == idx(2),1)-delay, d.spk, 'unif',0); 
  
  do_animate = false; 
  
  for pp = 1:length(d.spk) % identify original spike locations 

    if isempty(d.spk{pp}), continue, end

    t = d.spk{pp}(:,1) - delay; 
    x = d.spk{pp}(:,2);
    t_win = 0.1; 
    x_win = 8; 

    if do_animate, clf %#ok<UNRCH>
      rectangle('position',[t(1) 0 t_win 10],'FaceColor',[.8 .8 .8], ... 
                     'EdgeColor','none') 
      hold on, plot(t,x,'.', t(1),x(1),'o', t(1),x(1),'s')
      h = flipud(get(gca,'Children')); 
    end

    times = unique(t); 
    for tt = 1:length(times) % for each discrete timepoint 

      these = (t == times(tt));
      segs = cumsum(diff([-inf;x(these)])>2);
      prev = (t < times(tt) & t > times(tt)-t_win);

      if do_animate
        h(3).XData = t(these);  h(3).YData = x(these); %#ok<UNRCH>
        h(4).XData = t(prev);   h(4).YData = x(prev);
      end

      for ii = 1:max(segs)

        sel = these; 
        sel(these) = (segs == ii);
        adj = reshape(x(sel)+(-x_win:x_win),[],1); 

        if do_animate
          h(1).Position(1) = times(tt)-t_win; %#ok<UNRCH>
          h(1).Position(2) = min(adj);
          h(1).Position(4) = range(adj);
          pause(0.01);
        end

        adj = ismember(x(prev), adj);
        if ~any(adj)
          f.spikes.init{pp}(end+1,:) = [times(tt) round(median(x(sel)))+1];
        end
      end
    end
  end

  if pp == 1
    f.spikes.time = f.spikes.time{1};
    f.spikes.node = f.spikes.node{1};
    f.spikes.init = f.spikes.init{1};
    f.spiketimes  = f.spiketimes{1};
  end
else
  % determine spiketimes from Vm profile on "idx_prop" 
  f.spiketimes = tools.peakseek(f.V_example(:,2),0.1,-30); 
  f.spiketimes = f.time(f.spiketimes); 
end

save(cache_file,'-struct','f')

% a text file format was originally specified but has been removed, 
% .mat files are easy and the .dat format wasn't being maintained. 
% --- CDE 27-nov-20

return

function log_error(tid,iid,eid,vars,nrn_out)

named = evalin('caller','named');
if ~any(named('-log-e')), return, end

e_file = tools.file('sub~/axons/err-logs/');
if ~exist(e_file,'dir'), mkdir(e_file), end
e = fopen(sprintf('%serr_i%d_t%d_%s.log',e_file,iid,tid,eid),'wt'); 

fprintf(e,'ERROR LOG : %s\n',eid);
fprintf(e,'Thread # %d\n parfor index %d\n\n',tid,iid);

fprintf(e,'Output from NEURON: \n\n%s\n', nrn_out);

fprintf(e,'\nparameters { \n');
for ii = 1:size(vars,2)
    if ischar(vars{2,ii}), fprintf(e,'  %s ''%s''\n', vars{:,ii});
    elseif isstruct(vars{2,ii})
      f = fieldnames(vars{2,ii}); 
      f = sprintf(',%s',f{:});
         fprintf(e,'  %s struct [%s]\n',vars{1,ii},f(2:end));
    else fprintf(e,'  %s %.8f\n', vars{:,ii});
    end
end
fprintf(e,'}\n\nEND ERROR LOG\n');
fclose(e);


