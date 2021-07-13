
function [wave,wave_spk1] = spike_to_wave(g_id,time,sensors,g_xy,options)
% models.spike_to_wave retreives the saved membrane current file for a
% specified axon (g_id), and simulates the extracellular potential 
% 
% g_id = axon template index
% time = simulation time-points
% sensors = sensitivity (function handle)
% g_xy = xy co-ordinates of axons (units to match sensors) 
% % if g_xy is 2D, assume straight fascicles; 

% options
%   .raster = models.random_raster % what spike pattern to simulate?
% 
%   .dataroot / .class % what membrane current profiles to use? 
% 
%   .reference_z = 0 % what point in space are spikes aligned to?
%   .padding = 1     % left/right border padding on sensitivity function
%   .evoked_potential = [] % ? 
%   .efferent = 0 % flip propegation direction
%   .distortion = [1 1 1] % edit CV/Z-width/T-width
%   .get_all = 0 % return partial sum of the wave for each spike in _spk1
%   .plot_image = 0 % Make graphical output
%
% version 0.5 - 10 July 2020 Calvin Eiber
% version 0.4 - 27 May 2020  Calvin Eiber) 

if ~exist('options','var'), options = struct; end
if ~isfield(options,'reference_z'), options.reference_z = 0; end
if ~isfield(options,'padding'), options.padding = 1; end
if ~isfield(options,'evoked_potential'), options.evoked_potential = []; end
if ~isfield(options,'get_all'), options.get_all = false; end
if ~isfield(options,'plot_image'), options.plot_image = false; end
if isfield(options,'axons_folder'), options.dataroot = options.axons_folder; end

% There's something fishy happening in d.length, maybe an off-by-1 error?

nE =  numel(sensors); 
wave = zeros(length(time),nE);
wave_spk1 = zeros(length(time),nE,0);

FLAG_explicitSpikeTimePlace = false; 

if isfield(options,'raster')
  
  if isfield(options.raster,'spike') % explicit spike,node list supplied
    
    ax_g = (options.raster.axon_group == g_id); % Group mask    
    ax_t = ~cellfun(@isempty, {options.raster.spike(ax_g).time});
    if ~any(ax_t), return, end % escape early if no spikes in this group
    ax_g(ax_g) = ax_t; % apply existance mask to group mask
    
    if ismatrix(g_xy), g_xy = g_xy(ax_t,:); % fluff out g_xy
    else               g_xy = g_xy(ax_t,:,:); % fluff out g_xy [3D]
    end

    g_id = floor(g_id); % remove fascicle code (gg.ff) if relevent
    ax_n = (1:size(g_xy,1)); % fill in new indices
    
    options.raster.spike(~ax_g) = []; 
    options.raster.dt_dx(~ax_g) = []; 
    FLAG_explicitSpikeTimePlace = true; 
  
  elseif isfield(options.raster,'axon_group')

    ax_g = (options.raster.axon_group == g_id); % Group mask
    ax_k = ax_g(options.raster.spk_axon); % which times are in group?
    if ~any(ax_k), return, end % escape early if no spikes in this group
    
    ax_n = 0*options.raster.axon_group; % indexer
    ax_n(ax_g) = 1:sum(ax_g); % fill in new indices

    spike_time = options.raster.spk_time(ax_k); % get spiketimes      
    ax_n = ax_n(options.raster.spk_axon(ax_k)); % get new indices
    if ismatrix(g_xy), g_xy = g_xy(ax_n,:); % fluff out g_xy
    else               g_xy = g_xy(ax_n,:,:); % fluff out g_xy [3D]
    end

    g_id = floor(g_id); % remove fascicle code (gg.ff) if relevent
    clear ax_g ax_k
  else
    error simple_raster_input % (groups precomputed or no grouping)
  end
else 
  nrw = true; % do a warning
  if ~isfield(options,'jitter'), options.jitter = 0; nrw = false; end
  if ~isfield(options,'fraction'), options.fraction = 1; nrw = false; end
  if nrw, warning('ViNERS:Spike_to_wave:noRaster', ... 
                 ['No raster was supplied, please specify options.raster. ' ...
                  'Running spike test pattern using .fraction, .jitter'])
  end
  
  ok = (rand(size(g_xy,1),1) <= options.fraction);   
  if ismatrix(g_xy), g_xy = g_xy(ok,:); % fluff out g_xy
  else               g_xy = g_xy(ok,:,:); % fluff out g_xy [3D]
  end  
  
  ax_n = find(ok);   
  spike_time = randn(size(g_xy,1),1) * options.jitter;
  % if rand > options.fraction, continue, end
end

% For each spike, there are two quantities: 
%   spike_time (n x 1), g_xy (n x 2) or (n x 3 x nZ)

if isfield(options,'dataroot')
     d = sprintf('%s%s/n%03d_out.mat',options.dataroot,options.class,g_id);
elseif isfield(options,'class')
     d = tools.file(sprintf('axons~/%s/n%03d_out.mat', ... 
                                         options.class,g_id));
else d = sprintf('%sn%03d_out.mat',tools.cache('PATH'),g_id);
end

d = load(d); 

if isfield(options,'distortion')
  distortion_factor = options.distortion; 
  if isstruct(distortion_factor),     
      distortion_factor = [options.distortion.v ...
                           options.distortion.z ...
                           options.distortion.t];
  end
else distortion_factor = [1 1 1];
end
if isfield(options,'efferent') && options.efferent
   distortion_factor = [-1 -1 1] .* distortion_factor;
   % TODO: double-check this in a visual !!!!!
   % The idea is that the velocity and spatial patterns (but not the
   % temporal pattern) are inverted  
end
    
nZ = (numel(d.index)-1);
nP = (numel(d.length)-1) / nZ + 1;

t0 = find(d.V_example(:,2) > -20, 1); % from +models/analysis_interface.hoc AP_threshold = -20 (mV)
if isempty(t0), [~,t0] = max(d.V_example(:,2)); end
[~,p0] = min(abs(d.length(1:nP) - options.reference_z));

time = reshape(time,[],1);
% d.dt_dx is the delta-ms per length sample 

[Z,~] = sort(d.length);

i2v_list = cell(size(g_xy(:,1)));

dx = mean(diff(time)) / d.dt_dx; % integration scale for cumsum() 


for ii = 1:size(g_xy,1) % for each spike in raster ... 
  % printInfo('Group %d [u%03d]', gg,ii)
  
  c = d.I_membrane(:,2:end);
  % c(c<0) = c(c<0) / sum(c(c<0)) * -sum(c(c>0)); % Balance currents
  v = d.dt_dx;

  if FLAG_explicitSpikeTimePlace
    
    if isfinite(options.raster.dt_dx(ii)) && options.raster.dt_dx(ii) > 0      
      del_v = d.dt_dx ./ options.raster.dt_dx(ii);
      % del_v is < 1 if the observed spike is slower, and >1 if faster
      % we are assuming that slower spikes are also wider in time 
    else del_v = 1; 
    end

    if isempty(distortion_factor), tc = (d.time-d.time(t0));
      current = @(t) interp1(t+tc/del_v,c, time,'Linear',0); % at node n, current at each time
    else tc = (d.time-d.time(t0))/distortion_factor(3);
      current = @(t) interp1(t+tc/del_v,c, time,'Linear',0);
    end
    
  else
    
    if isempty(distortion_factor), t = spike_time(ii) + (d.time-d.time(t0));
      current = @(n) interp1(t, c, time-n*v,'Linear',0); % at node n, current at each time
    else t = spike_time(ii) + (d.time-d.time(t0))/distortion_factor(3);
      current = @(n) interp1(t,c,time-(n*v/distortion_factor(1)),'Linear',0);  
    end
    if ~isempty(options.evoked_potential), 
      poststimulus = 1./(1+exp((spike_time(ii)-time+options.evoked_potential) ...
                                  / options.raster_opts.tau1 * 3));
      current = @(n) current(n) .* poststimulus;
    end

    clear c v t
  end
  
  %%
  
  %%  
  if options.plot_image
    %% DEBUG image to check prop velocity and spike timing
    % this is basically a panel of the figure generated by
    % distortion_visualisation
    distortion_image = [];
    ee = 2; 
    if any(ax_n(1:ii-1) == ax_n(ii)),
         i2v = i2v_list{find(ax_n(1:ii-1) == ax_n(ii),1)}(:,ee);
    else
      if isempty(distortion_factor), i2v = get_i2v(sensors,ee,Z,g_xy(ii,:,:));
      else i2v = get_i2v(sensors,ee,Z/distortion_factor(2),g_xy(ii,:,:));
      end
    end
    % tx = tic;
    for pp = 1:nP-1
      i2v_seg = i2v( (pp-1)*nZ+(1:nZ) );
      distortion_image = [distortion_image current(pp-p0) * i2v_seg];
    end
    % tx = toc(tx);
  
    if ~exist('spike_time','var'), spike_time(ii) = 0; end
      
    
    clf, % subplot(3,1,[1 2])
    imagesc(time, Z, log10(abs(distortion_image))')
    hold on, plot(spike_time(ii),options.reference_z,'r+','markersize',30)
    xlabel('time, ms'), ylabel('Z, mm')    
    
    plot(time,sum(distortion_image,2),'-w','LineWidth',1.1)
    % error('Generated example distortion image')
    return
  end
  
  %%
  if any(ax_n(1:ii-1) == ax_n(ii)),
    i2v = i2v_list{find(ax_n(1:ii-1) == ax_n(ii),1)};
  else
    for ee = 1:nE % for each electrode
      if isempty(distortion_factor), i2v(:,ee) = get_i2v(sensors,ee,Z,g_xy(ii,:,:));
      else i2v(:,ee) = get_i2v(sensors,ee,Z/distortion_factor(2),g_xy(ii,:,:));
      end
    end    
    i2v(isnan(i2v)) = 0; 
    
    % options.i2v_taper = 2*nZ; 
    % gtaper = sin(linspace(0,pi/2,options.i2v_taper)');
    % i2v(1:options.i2v_taper,:) = i2v(1:options.i2v_taper,:) .* gtaper;
    % i2v(end:-1:end-options.i2v_taper+1,:) = i2v(end:-1:end-options.i2v_taper+1,:) .* gtaper;
    
    i2v_list{ii} = i2v;
  end

  if mean(i2v(:) == 0) > 0.5, continue, end    
  
  %%
  
  V = 0*wave; 
  balance_ = @(x) -x(x<0).*sum(x(x>0))./sum(x(x<0)); 
  
  if FLAG_explicitSpikeTimePlace
    
    spk = options.raster.spike(ii);
    
    for tt = 1:numel(spk.time)
      
      pp = spk.node(tt)+1;
      
      if tt > 1
        pre = abs(spk.node(1:tt-1) - spk.node(tt)) < 3;
        pre = spk.node(tt) - spk.node(find(pre,1,'last'));
        dir_fwd = any(pre > 0);
        dir_rev = any(pre < 0);
      else
        dir_fwd = true;
        dir_rev = true;
      end
      
      C0 = current(spk.time(tt)); 
      
      if dir_fwd
        idx = (pp-1)*nZ+(1:nZ); 
        ok  = (idx <= numel(Z)) & idx > 0;
        i2v_seg = i2v( idx(ok), : );
        V = V + C0(:,ok) * i2v_seg;
      end
      if dir_rev
        idx = (pp-1)*nZ-(1:nZ)+2;
        ok  = (idx <= numel(Z)) & idx > 0;
        i2v_seg = i2v( idx(ok), : );
        V = V + C0(:,ok) * i2v_seg;
      end
    end
    
    wave = wave + V ;

    
    if ii == 1,  wave_spk1 = wave;
    elseif options.get_all % aka options.extract_all_waves
      wave_spk1(:,:,end+1) = ( V ); %#ok<AGROW>
    end
    
    continue
    
  end % explicit spike-time,location tuples (e.g. ECAP simulation)
  
  
  %%
  if ~isempty(options.padding)

    C0 = current(0);  C0(C0<0) = balance_(C0);
    i2v_seg = ones(nZ,1)*i2v(1,:);
    
    edge_effect = cumsum(C0 * i2v_seg,'reverse');
    adi = mean(abs(diff(edge_effect([1 1:end],:))),2);     
    aoi = (adi > 1.5 * median(adi(adi>0)));
    nPad = ceil(sum(aoi) * dx * options.padding);
    vPad = 0*V;

    for pp = -nPad:0, % LHS padding
      C = current(pp-p0);
      % C(C<0) = balance_(C); 
      vPad = vPad + (C * i2v_seg);
    end
    
    i2v_seg = ones(nZ,1) * i2v(end,:);
    for pp = nP+(1:nPad+1), % RHS padding
      C = current(pp-p0);
      % C(C<0) = balance_(C); 
      vPad = vPad + (C * i2v_seg);
    end
    
    V = V + vPad;

  else nPad = 0; vPad = 0; 
  end

  
  %% Spatial summation of fields as the spike passes over the electrode
  
  for pp = 1:nP-1 % insert points
      C = current(pp-p0);
      C(C<0) = balance_(C); 
      i2v_seg = i2v((pp-1)*nZ+(1:nZ),:);
      V = V + C * i2v_seg;
  end

  % last segment - extrapolate i2v
  i2v_seg = i2v((pp-1)*nZ+(1:nZ),:) - i2v(pp*nZ-nZ+1,:) + i2v(pp*nZ+1,:);
  V = V + current(pp-p0+1) * i2v_seg;

  %% Apply some corrections to get rid of edge artifacts 
  % (not sure how effective this lot actually is)
  
  % C0 = current(0-nPad-p0-1.5);  % C0(C0<0) = balance_(C0);
  % C1 = current(nP+nPad+2.5-p0); % C1(C1<0) = balance_(C1);
  % edge_effect = cumsum(C1 * V2I_seg) + cumsum(C0 * i2v(1:nZ),'reverse'); 
  % edge_effect = edge_effect - linspace(edge_effect(1), ...
  %               edge_effect(end),numel(edge_effect))';

  taper = max(0, abs(time-spike_time(ii)) - (d.dt_dx * nP));
  taper = exp(-taper.^2);

  roi = (taper > 1e-2 & taper < 0.5);
  % trend = fit(time(roi),V(roi),'poly3');
  trend = (time(roi) * [1 0] + [0 1]) \ V(roi); 
  roi = (taper > 1e-3); 
  % trend = trend(time(roi));
  trend = time(roi) * trend(1) + trend(2);
  V(roi) = V(roi) - trend;

  if any(isnan(V(:)))    
    error('NaN in computed V_extracellular')
  end
    
  wave = wave + ( V .* taper );
% wave(:,ee) = wave(:,ee) + (V+edge_effect/k) .* taper;

  if isfield(options,'plot_wave') && ii == options.plot_wave 
    %%
    clf, hold on, lc = lines(7);
    plot(time,V(:,end)-vPad(:,end),'--','color',(lc(1,:)+1)/2)    
    plot(time,V(:,end),'-','color',lc(1,:),'linewidth',1)

    plot(time,vPad(:,end),'-','Color',lc(6,:))

    if exist('C1','var')      
      plot(time,-edge_effect(:,end)/dx,'-','Color',(lc(2,:)+2)/3,'linewidth',1)
      plot(time,(V+edge_effect(:,end)/dx) .* taper,'k-','LineWidth',1.2)
    else plot(time,V(:,end) .* taper,'k-','LineWidth',1.2)        
      plot(time(roi),trend(:,end),'-','Color',(lc(2,:)+2)/3,'linewidth',1)
    end

    plot(time,0.15 * taper,':','LineWidth',1.3,'Color',lc(5,:))
    xlim(spike_time(ii) + [-10 10])
  end

  if ii == 1,  wave_spk1 = wave;
  elseif options.get_all % aka options.extract_all_waves
      wave_spk1(:,:,end+1) = ( V .* taper ); %#ok<AGROW>
  end
end % loop over units

% WAVE = µV = (nA) / (4*pi*sigma) (ohm.m) / distance (mm) 
% the mm and the nA cancel to yield µV, no conversion necessary

function i2v = get_i2v(sensors, ee, Z, xy)

if ndims(xy) == 3
 
  xy = permute(xy,[3 2 1]);
  u = sqrt(sum(diff(xy).^2,2)); 
  u = cumsum(conv(u([1 1:end end]),[1;1]/2,'valid'));
  [~,z0] = min(abs(xy(:,1))); u = u-u(z0);
  
  Zq = fliplr(interp1(u,xy,Z));
  i2v = sensors{ee}(Zq(:,1),Zq(:,2),Zq(:,3));
  
else 
  
  i2v = sensors{ee}(0*Z+xy(1), 0*Z+xy(2), Z);
  if all(isnan(i2v(Z >= 0))) && ~all(isnan(i2v(Z <= 0)))    
    % To save time and space, only the right half (+z) was simulated.
    %   Use symmetry to get the other half of the sensitivity curve. 
    i2v(Z<0) = sensors{end-ee+1}(0*Z(Z<0)+xy(1), 0*Z(Z<0)+xy(2), -Z(Z<0));
  end
end

ok = ~isnan(i2v); % Patch over discontinuities
if mean(ok) < 0.5, i2v(:) = 0; return, end     
if all(ok), return, end
i2v = interp1(Z(ok),i2v(ok),Z,'linear','extrap');

% Construct tapered i2v function [not wanted!]
% i2v_edge = max(0, i2v(1) - (nPad:-1:0)*(i2v(1,2)-i2v(1)));
% i2v_edge = i2v_edge .* (cos(linspace(pi/2,0,nPad+1)).^2);     
% i2v_edge = max(0, i2v(1,end) - (0:nPad) * (i2v(1,end-1)-i2v(1,end)));
% i2v_edge = i2v_edge .* (cos(linspace(0,pi/2,nPad+1)).^2); 



