function export_NSx_file(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

%% Load an example NS5 file recorded at the Bionics Institute to emulate

if isempty(which('openNEV')) % BlackRock Matlab API
  p = tools.configuration('NPMK'); % Ask pelvic-nerve-model nicely
  if isempty(p) && ~isempty(which('expdir')), p = expdir('NPMK'); end
  if isempty(p), p = fullfile(userpath,'Add-Ons','NPMK'); end
  path(path,p)
  path(path,[p filesep 'NSx Utilities'])
  clear p
end

%%

opts = set_output_configuration(varargin{:}); 

flag_UNIT_TEST = 0; 

if flag_UNIT_TEST
  
  disp('Unit-Test') %#ok<UNRCH>
  
  opts.input = {tools.file('~\primary\sub-2\waves\long runs\epoch_k1_f2.9_c0.3 (1).mat')};

  opts.output(2) = opts.output(1);
  opts.output(2).seg.input_noise = 30;  
  
  opts.output(2).seg(2) = opts.output(2).seg(1); 
  opts.output(2).seg(2).input_noise = 80; 
end

indices.file  = [];
indices.f_out = [];
indices.s_out = [];
indices.p_out = []; 

for index = 1:numel(opts.output)
  for ss = 1:numel(opts.output(index).seg)
    
    fidx = reshape(opts.output(index).seg(ss).fileIndex,1,[]); 
    indices.file = [indices.file fidx];
    indices.f_out = [indices.f_out 0*fidx + index];
    indices.s_out = [indices.s_out 0*fidx + ss];
    indices.p_out = [indices.p_out (1:numel(fidx))];
  end
end


defaults.channels = 'default';
if any(named('-pair')), defaults.channels = get_('-pair'); end
if any(named('-elec')), defaults.channels = get_('-elec'); end

%%
waves = cell(size(indices.file)); 

for index = 1:numel(opts.input)

  if ~any(indices.file == index), continue, end
  
  if any(opts.input{index} == '~'), 
         opts.input{index} = tools.file(opts.input{index}); 
  end
  
  if ~exist(opts.input{index},'file'), 
    error('Could not find "%s"', opts.input{index})
  end
  
  SIM = load(opts.input{index});
  SIM.fs = 1./mean(diff(SIM.time)); 

  if abs(SIM.fs - 30) > 1e-10
    warning('ViNERS:exportNSx:resampleTime', ...
            'The imported data has a Fs of %0.1f kHz, resampling to 30k', SIM.fs)
    error TODO_implement
  end

  nE = size(SIM.waves,2); 
  
  if ~any(named('-no-detrend'))
    roi = abs(SIM.time) < max(abs(SIM.time)) - 25;
    SIM.waves = tools.detrend_wave(SIM.waves,SIM.time,roi);
  end
  
  for oidx = find(indices.file == index)
    
    ff = indices.f_out(oidx);
    ss = indices.s_out(oidx);
    pp = indices.p_out(oidx);
    
    if ischar(opts.output(ff).channels)
      if ischar(defaults.channels), 
        opts.output(ff).channels = reshape( 1:nE,2,[] )';
      else opts.output(ff).channels = defaults.channels;
      end
    end, chan = opts.output(ff).channels; 

    if max(chan(:)) > nE
      error('electrode %d requested, but only %d electrodes found in "%s"', ...
            max(chan(:)), nE, opts.input{index})
    end
    
    roi = opts.output(ff).seg(ss).input_roi; 
    if numel(roi) > 1, roi = roi(pp); end
    if roi < 0, roi = max(abs(SIM.time)) + roi; end
    roi = abs(SIM.time) < roi; 
    
    delt = opts.output(ff).seg(ss).input_shift; 
    if numel(delt) > 1, delt = delt(pp); end
    if delt ~= 0
      delt = round(delt * SIM.fs);
      roi = circshift(roi,[1 -delt]);
    end
    
    aidx = opts.output(ff).seg(ss).axonClass(pp); 
    
    for cc = 1:size(chan,1)

      if ndims(SIM.waves) == 3
        wave = SIM.waves(roi,chan(cc,:),aidx); % fascicles pre-summed
      elseif ischar(opts.output(ff).seg(ss).fascicles)
        wave = sum(SIM.waves(roi,chan(cc,:),:,aidx),3);
      else sel = opts.output(ff).seg(ss).fascicles{pp};
        sel(sel > size(SIM.waves,3)) = []; 
        wave = sum(SIM.waves(roi,chan(cc,:),sel,aidx),3);
      end

      if size(wave,2) > 1
        wave(:,2:end) = - wave(:,2:end) ./ numel(wave(1,2:end)); 
        wave = sum(wave,2);
      end
      
      n_uV = opts.output(ff).seg(ss).input_noise; 
      if numel(n_uV) > 1, n_uV = n_uV(pp); end
    
      gain = opts.output(ff).seg(ss).input_noise; 
      if numel(gain) > 1, gain = gain(pp); end
    
      wave = wave * gain + randn(size(wave)) * n_uV; 
      waves{oidx}(:,cc) = wave; 
    end
  
    sel = (indices.f_out == ff);    
    
    if any(cellfun(@isempty,waves(sel))), continue, end
    
    write_NEV_file(opts.output(ff),waves(sel),varargin{:});
    [waves{sel}] = deal('done'); 

  end
  
end

% Data processing options:
% -noise 10 (µV) : define Gaussian noise amplitude
% -roi 25 (ms)   : define time window to exclude edge effects 
%                   (symmetric around 0 or 2-elementvector)
% -elec [...]    : define recording electrode pairs 
%                   (synonymous with -chan, -pair)
% -fasc [1:nF]   : reduce recording to the specified fascicle IDs

% assert(  1./mean(diff(D.time)) == 30 )


%% Default output configuration
function s = set_output_configuration(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};
f_ = @(f) [f.folder filesep f.name];

if any(named('-co')), s = get_('-co'); else s = struct; end 
% -config or -composite

if ~isfield(s,'input') % Parse inputFile from arg{1} or select
  
  s.input = {};
  
  if nargin > 0 && ischar(varargin{1}) 
    if any(varargin{1} == '~'), varargin{1} = tools.file(varargin{1});
    end
    if ~exist(varargin{1},'file'), varargin{1} = dir(varargin{1}); end
  end
  if any(named('-in')),  s.input = get_('-in');
  elseif nargin > 0 && (~ischar(varargin(1)) || exist(varargin{1},'file'))
                         s.input = varargin{1};
  end
  if isstruct(s.input) && isfield(s.input,'folder')
    s.input = arrayfun(f_,s.input,'unif',0);
  end
  if ~iscell(s.input), s.input = {s.input}; end
  
  if isempty(s.input), 
    s.input = {tools.file('sub~/waves/*.mat','-prompt')};
    
    if any(named('-all')) || any(named('-get'))
      [fp,~] = fileparts(s.input{1}); 
      s.input = arrayfun(f_,dir([fp '/*.mat']),'unif',0); 
    end
  end
end

if isfield(s,'output'), return, end

s.output = []; 

% The default composition is to copy the input mat files to the output nsx
% files with no substantive changes. gather optional user input to set up 
% default values (in "d"):  

d.noise = 10;
d.time  = -25;
d.shift = 0;
d.gain  = 1; 

d.axonClass = 1:4; 
d.fascicle = 'all'; 

if any(named('-ax')),    d.axonClass = get_('-ax'); end
if any(named('-noise')), d.noise = get_('-noise');  end
if any(named('-roi')),   d.time = get_('-roi');     end
if any(named('-shift')), d.shift = get_('-shift');  end
if any(named('-gain')),  d.gain = get_('-gain');    end
 
size_check_ = @(u) all(size(u) == 1) || ...
                   all(size(u) == size(d.axonClass)); 
for prop = {'noise','time','gain','shift'}
  if size_check_(d.(prop{1})), continue, end
  error('%s was an invalid size (must be 1x1 or %dx%d, was %dx%d)', ... 
            prop{1},size(d.axonClass), size(d.(prop{1})))
end


% model.waves is a nTime x nChannels x nFascicles x nTypes 
%             OR a nTime x nChannels x nTypes matrix in units of uV

% s.output.seg(nV).fileIndex = [1 x nW] index into s.input
% s.output.seg(nV).fascicles = 'all' or {1 x nW} list ? 
% s.output.seg(nV).axonClass = [1 x nW] 
% s.output.seg(nV).input.shift = scalar [0] or [1 x nW] in // ms //
% s.output.seg(nV).input.gain  = scalar [1] or [1 x nW]
% s.output.seg(nV).input.noise = scalar [10] or [1 x nW] in // uV //
% s.output.seg(nV).input.roi   = scalar [-25] or [1 x nW] in // ms //
% s.file.fileName = 'whatever.ns5'

for ff = 1:numel(s.input)

  this = struct; 
  this.filename = sprintf('simulated_data_%02d_%%d.*',ff);
  this.seg.fileIndex = ff * ones(size(d.axonClass));
  this.seg.fascicles = d.fascicle;
  this.seg.axonClass = d.axonClass;
  this.seg.input_shift = d.shift; 
  this.seg.input_gain  = d.gain;
  this.seg.input_noise = d.noise;
  this.seg.input_roi   = d.time; 
  
  this.channels = 'default';
  
  if isempty(s.output), s.output = this;
  else s.output(end+1) = this;
  end  
end

return

%% Convert waves to NEV file and write it
function write_NEV_file(this,waves,varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

persistent SRC
if isempty(SRC)

  if any(named('-nev-')), f_example = get_('-nev-');
  elseif any(named('-refer')), f_example = get_('-refer');
  else f_example = [expdir('Rat12_016') 'Baseline005.ns5']; 
  end, SRC = openNSx(f_example); 

end

NS5 = SRC;
NS5.Data(:) = []; 

widx = arrayfun(@(s) s*ones(size(this.seg(s).fileIndex)), ...
                    1:numel(this.seg),'unif',0);
widx = [widx{:}];
assert(numel(widx) == numel(waves))

%%
nC = max(cellfun(@(w) size(w,2), waves));
nE = numel(NS5.ElectrodesInfo);

if nE > nC && any(named('-no-fill'))
  warning('ViNERS:exportNSx:addingChannels', ...
            '%s %d channels (template had %d)', ...
            'The exported NS5 file will contain', nC, nE)
    error TODO_implement        
elseif nC > nE
  warning('ViNERS:exportNSx:addingChannels', ...
            '%s %d channels (template had %d)', ...
            'The exported NS5 file will contain', nC, nE)
    error TODO_implement
end

a2d = zeros(nE,2);

for cc = 1:nE
  
  val = double([NS5.ElectrodesInfo(cc).MinDigiValue ...
                NS5.ElectrodesInfo(cc).MaxDigiValue ...
                NS5.ElectrodesInfo(cc).MinAnalogValue ...
                NS5.ElectrodesInfo(cc).MaxAnalogValue]);      
  a2d(cc,:) = [(val(4)+val(3))/2 range(val(1:2))/range(val(3:4))];
end

%% Gather waves for each segment

for ss = 1:numel(this.seg)
  
  sel = widx == ss; 
  min_nT = min(cellfun(@(w) size(w,1), waves(sel)));
  
  for pp = find(sel)  
    waves{pp} = waves{pp}(1:min_nT,:); 
  end
  
  wave = sum(cat(3,waves{sel}),3);
  
  for cc = 1:nE
    
    if cc > nC
      wave(:,cc) = randn(size(wave(:,1))) * this.seg(ss).input_noise; 
    end
    
    wave(:,cc) = (wave(:,cc) + a2d(cc,1)) * a2d(cc,2); 
  end
  
  waves = apply_NSx_filters(NS5,waves); 
  
  NS5.Data = [NS5.Data cast(wave','like',NS5.Data)];
end

%% Set up metadata correctly

NS5.MetaTags.DateTime = datestr(now);
NS5.MetaTags.DateTimeRaw = cast([round(datevec(now)) 0 0],'like',...
                                    NS5.MetaTags.DateTimeRaw); 

NS5.MetaTags.DataPoints = size(NS5.Data,2);
NS5.MetaTags.DataDurationSec = size(NS5.Data,2) / 3e4;
NS5.MetaTags.DataPointsSec = size(NS5.Data,2) / 3e4;

comment = 'simulated_data from pelvic-nerve-model eiber et al 2020';
comment(end+1 : length(NS5.MetaTags.Comment)) = 0; 
NS5.MetaTags.Comment = comment;

this.filename = tools.file('get',strrep(this.filename,'.*','.ns5'),'next');

if isempty(fileparts(this.filename)),
  this.filename = ['.' filesep this.filename];
end

saveNSx(NS5,this.filename,'-y')

return

%% Apply NS5 filter settings to waves

function waves = apply_NSx_filters(this,waves)


if iscell(waves), 
  waves = cellfun(@(w) apply_NSx_filters(this,w), waves, 'unif',0); 
  return
end

if this.MetaTags(1).HighFreqCorner == 0, return, end % No filter set

error TODO_filtering

return