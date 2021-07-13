
function raster = random_raster(varargin)
% Generate a random raster of axon spikes following the specified
% non-stationary population firing rate. Two syntaxes are supported: 
% 
% raster = random_raster(time, g_xy, [opts]) 
% raster = models.random_raster( 'name',value, ... ) [NOT IMPLEMENTED YET]
% 
% Several non-stationary population codes are supported (as detailed below)
%  as well as statonary (poisson) spiking. Within a population of the
%  specified size, individual axons' spike-rates are lognormally
%  distributed around the population average spike-rate. 
%
% The following population codes are supported: 
%  .type = 'poisson' (or 'stationary' or 'flat'): rate(t) = base
%  .fb [2]  = base : baseline, imp/s (constant rate)
%
%  .type = 'raised_cosine (cos)' : rate(t) = 
%  (peak - base) * ( cos(2 pi freq t + phi)/2 + 1/2 )^expo + base
% 
%  .fb [2]  = base : baseline, imp/s
%  .fp [40] = peak : max spikerate, imp/s 
%  .ph [0]  = phi  : phase offset, radians
%  .fc [30] = freq : frequency of population firing-rate modulation, Hz
%  .ex [4]  = expo : exponent 
% 
% .type = 'double_exponential (exp2)' : rate(t) = base, t < t0 ...
%  (peak-base) * (1-e^( (t0-t)/tau_1 )) * e^( (t0-t)/tau_2 ) + base, t > t0
% 
%  .fb [2]  = peak : baseline, imp/s
%  .fp [40] = base : peak, imp/s 
%  .t0 [0]  = t0   : onset time, ms
%  .tau1 [] = tau_1 : onset time constant, ms
%  .tau2 [] = tau_2 : decay time constant, ms
% 
% .type = 'step_function (step)' : rate(t) = 
%    (end - start) / ( 1 + e^( (t0 - t) / tau ) ) + start
% 
%  .fb [2]  = start : starting rate, imp/s
%  .fp [40] = end   : ending rate, imp/s 
%  .t0 [0]  = t0    : half-maximal time, ms
%  .ex [4]  = tau   : 1 / slope of step at half-maximal point, ms
%                     (the step covers ~46% of the difference between peak 
%                      and base in from t0-tau to t0+tau)
% Other options: 
% 
%  .dt [1] bin size for spike-rate histogram
%  .ax_sd [1] imp/s, stdev * mean rate increasing  
% 
% Othe options (name,value only): 
% 
% -mat  : save output as .mat file
% -xml  : save output as .xml file
% -json : save output as .json file
% -roi [start end] or [0 end] (name-value only)
% -count : number of axons in population
% 
% v0.3 CDE 16 June 2021

if nargin > 1 && isnumeric(varargin{1})
  raster = make_random_raster(varargin{:}); 
  return
end

varargin = tools.opts_to_args(varargin,'raster');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

% I haven't decided how to implement any of this yet. 

error TODO_generate_and_save_raster_XML











%%

function raster = make_random_raster(time, g_xy, opts)

% set up default values for opts: 
default.ty = 'raised_cosine';
default.fb = 2; % imp/s, baseline
default.fp = 40; % imp/s, peak
default.fc = 30; % Hz, oscillation
default.ex = 4; % exponent
default.ph = 0; % phase, radians
default.pop_dist = 'lognormal'; % poisson / exponental? 
default.ax_sd = 1; % imp/s, stdev * mean rate
default.dt = 1; % ms bin size for rate 

default.tau1 = 1; % for 'double_exponential'
default.tau2 = 10;

if nargin < 3, opts = struct; end

for f = fieldnames(default)'
  if ~isfield(opts,f{1}), opts.(f{1}) = default.(f{1}); end
end

if isfield(opts,'type'), opts.ty = opts.type; end
if isfield(opts,'t0'), opts.ph = opts.t0; end

%%
nA = size(g_xy,1); 

if isa(opts.ty,'function_handle'), response_fcn = opts.ty;
else
 switch opts.ty
  case {'flat','poisson','stationary'}
      response_fcn = @(t,p) ones(size(t))*p.fb;
  case {'raised_cosine', 'cos'}
    response_fcn = @(t,p) (p.fp - p.fb)*((cos(2*pi*t*p.fc/1000 + ...
                                      p.ph)+1)/2).^p.ex + p.fb;
  case {'double_exponential','exp2'}
    response_fcn = @(t,p) (p.fp - p.fb) * (t > p.ph) .* ...
                             (1-exp(-(t-p.ph)./p.tau1)) .* ...
                               (exp(-(t-p.ph)./p.tau2)) + p.fb;
  case {'step_function','step','sigmoid'}
    response_fcn = @(t,p) (p.fp - p.fb) ./ ...
                             (1+exp(-(t-p.ph)./p.ex)) + p.fb;
                             
  otherwise error('Unknown response function %s', opts.ty)
 end
end

pop_SR = response_fcn(time,opts);

if isfield(opts,'modulation')
  mask = abs(opts.modulation(time));
  pop_SR = pop_SR .* mask + opts.fb .* (1-mask); 
end


nK = round( nA * mean(pop_SR) * range(time) / 1000); % total number of spikes

if strncmp(opts.pop_dist,'norm',4)
  ax_K = mean(pop_SR) + opts.ax_sd * mean(pop_SR) * randn(nA,1); % spikes per each axon
else
  ax_K = mean(pop_SR) * exp( opts.ax_sd .* randn(nA,1) ); 
end

ax_K(ax_K < 0) = 0; ax_N = round(ax_K ./ sum(ax_K) * nK); % clamp and renomalise

pop_csr = cumsum(pop_SR) / sum(pop_SR);
spike_time = rand(nK,1);
spike_axon = 0*spike_time; 

for kk = 1:nK

    spike_time(kk) = time(find(pop_csr >= spike_time(kk),1));
    
    if any(ax_N)
        sel = find(cumsum(ax_N)/sum(ax_N) > rand,1); 
        spike_axon(kk) = sel; 
        ax_N(sel) = ax_N(sel)-1; 
    else
        spike_axon(kk) = ceil(rand*numel(ax_N));
    end
end


bin_x = min(time):(opts.dt):max(time);
bin_y = hist(spike_time,bin_x); 
bin_y = bin_y ./ mean(diff(bin_x)) / nA * 1000;



raster.spk_time = spike_time;
raster.spk_axon = spike_axon;
raster.spk_rate = bin_y; 
raster.bin_time = bin_x; 
raster.bin_rate = response_fcn(bin_x,opts);
raster.pop_rate = pop_SR; 

if nargout > 0 && (~isfield(opts,'do_plot') || opts.do_plot), return, end % skip visualisation

clf, hold on, C = lines(7); 
bar(bin_x,bin_y,1,'FaceColor',[.7 .7 .7],'EdgeColor','none')
y_ = @(n) opts.fp * (n/max(n) + 1.02);
plot(spike_time,y_(spike_axon),'k.');
plot(time, pop_SR,'-s','LineWidth',1.1,'Color',C(1,:))
% plot(raster.bin_time, raster.bin_rate,'-s','LineWidth',1.1,'Color',C(1,:))
tools.tidyPlot


%%



