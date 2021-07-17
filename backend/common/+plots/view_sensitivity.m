
function view_sensitivity(varargin)
% Plots.view_sensitivity( filename, ... )
% 
% Produce a visualisation of fascicle recording sensitivity, stimulus 
%   v_extracellular, or activating function. The visualisation selected
%   depends on the number of fascicles and the number of electrodes. 
% 
% if # fascicles > 1 and # electrodes > 1: summary cross-sections 
%  (1 / electrode) for peak field, (optionally) field spatial width
% if # fascicles = 1 and # electrode > 1: longitudinal traces, clickable
%  (1 / electrode) cross-sections switch electrodes with a click on axes.
% if # fascicles = 1 and # electrodes = 1: clickable explorer view with
%  clickable section (no heatmap), longitidudinal trace, and selected point
%  extuded through the mesh showing the interpolation. 
% 
% additional figures can also be generated using flags: 
%   -geo:  Produce array geometry figure,
%   -xy:   Generate a scatterplot comparison from the summary
%           cross-section view. 
% 
% If called with no filename, select an eidors sensitivity or stimulus
%   file using a file chooser dialog.
% 
% Additional options: 
%  -bi:  use bipolar electrode montage for recording sensitivity
%  -af:  display computed activating function (Rattay 1999)
%  -mf:  override logic to make multi-fascicle section figure
%  -hw:  show spatial widths in multi-fascicle figure
%  -f:   set fascicles to display.
%  -e:   set electrodes to display.
%  -q:   quiet mode, less console output 
% 
% Version 0.6 13-Oct-2020 CDE

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

% files or gather from input
if nargin > 0, file = varargin{1}; else file = {}; end
if isempty(file) || file(1) == '-'
   file = {'*.mat'}; 
   if any(named('-af')), file = file([2 1]); end
  [file,folder] = uigetfile(file,[],tools.file('out~\'));
  if all(folder == 0), return, end % cancelled
  disp(file)
  file = strcat(folder,file);
elseif strcmpi(file,'auto')
end
if contains(file,'~'), file = tools.file(file); end

%%

EM = load(file);
EM.filename = file;

FLAG_sfap = isfield(EM,'axons') && isfield(EM,'eidors_file'); 
if FLAG_sfap, plots.view_sfap(varargin{:}), return, end

FLAG_stimulus = isfield(EM,'v_extracellular');
FLAG_multimesh = (numel(EM.model) > 1); 
FLAG_verbose = ~any(named('-q')) || isdeployed; 

%%
if FLAG_stimulus
  if FLAG_verbose
    if any(named('-AF'))
         disp('Formatting V_e into activating function in fascicles (-AF)')
    else disp('Formatting V_e into fascicles (-STIM)')
    end
  end
  EM = convert_stim2fascicles(EM);
elseif any(named('-bi'))
  if FLAG_verbose 
    disp('Formatting response into bipolar measurements')
  end
  EM = convert_mono2bipolar(EM); 
  
end

mk_shortcuts(EM,FLAG_multimesh)

%%

if any(named('-align-EM-data'))
    %% Align EM data to artificial geometry
    
    error TODO_import_EM_axons
    
    xy = [z_(1) y_(1)]; 
    f_outline = convhull(xy(:,1),xy(:,2));
    f_outline = xy(f_outline([1:end 1]),:);
    
    c_xy = mean([T.Centroid_X T.Centroid_Y]);
    xy = ([T.Centroid_X T.Centroid_Y] - c_xy)/1000;
    % A is roughly in descending A.Centroid_X order
    
    rotate = @(v, a) v * [cos(a) -sin(a); sin(a) cos(a)]; 
    xy = rotate(xy, pi/3) + mean([z_(fac_(f)) y_(fac_(f))]);

    clf, hold on
    plot(f_outline(:,1),f_outline(:,2),'-','Color',C(1,:),'LineWidth',1.2)
    plot(z_(1), y_(1), '.','Color',W(1,1.5))
    plot(xy(:,1),xy(:,2),'.','Color',C(2,:))
    axis equal tight, tools.tidyPlot
end

clear rotate c_xy f ee fn fp 

%% Generate requested plots 

if any(named('-geo')), array_geometry_figure(EM), return, end

%%

nE = size(EM.Fascicle1.pot,1); 
electrode_list = 1:nE;    % [1 2 9 10]; 

if any(named('-c')), electrode_list = get_('-c'); end
if any(named('-e')), electrode_list = get_('-e'); end

if any(electrode_list > nE) 
  warning('File x only has Y electrodes')
  electrode_list(electrode_list > nE) = []; 
end

opts.fascicle = 1:nF;
if any(named('-f')), opts.fascicle = get_('-f'); end
opts.electrode = electrode_list; 

opts.do_hwidth = any(named('-hw'));
opts.do_parametric = any(named('-par')); 
opts.do_image = any(named('-im')); 

opts.force_mf = any(named('-mf'));
opts.xy_res = 21; 
if any(named('-res')), opts.xy_res = get_('-res'); end
opts.pk_quantile = 1;
if any(named('-pk-q')), opts.pk_quantile = get_('-pk-q'); end

if any(named('-pca:')), opts.do_PCA = get_('-pca:'); % for -single-f
elseif any(named('-pca')), opts.do_PCA = 3; 
end

opts.debug_thinlayer_idx = any(named('-no-peri')); 

clf, choose_sensitivity_figure(EM,opts)

%%

if any(named('-xy')) && numel(electrode_list) > 1
%%
  
  idx = [1 2];  
  if any(named('-xy-i')), idx = get_('-xy-i'); end

  ax = flipud(get(gcf,'Children'));
  if any(arrayfun(@(x)strcmp(x.Tag,'Colorbar'),ax))  
    ax = ax(1:4:end); % Get just "peak" panels
  end


  x = [ax(idx(1)).Children(1:2:end).CData];
  y = [ax(idx(2)).Children(1:2:end).CData];

  figure, clf, plot(x,y,'o')
  tools.tidyPlot, axis equal square, hold on
  plot(xlim,xlim,'Color',[0 0 0 0.3])

  xlabel([ax(idx(1)).YLabel.String ' Peak (�V/�A)'])
  ylabel([ax(idx(2)).YLabel.String ' Peak (�V/�A)'])
  
end

%% Top-down view of array
function array_geometry_figure(EM) 

  get_ = evalin('caller','get_');
  named = evalin('caller','named');
  FLAG_multimesh = (numel(EM.model) > 1); 
  
  
  mk_shortcuts(EM,FLAG_multimesh);
  
  elec_sensor = cell(nE,1); 

  f_id = 1; 
  if any(named('-f')), f_id = get_('-f'); end
  fascicleN = sprintf('Fascicle%d',f_id);

  for ee = 1:nE
    elec_sensor{ee} = scatteredInterpolant(z_(f_id), y_(f_id), x_(f_id), ...
                                   EM.(fascicleN).pot(ee,:)','linear','none'); 
  end
  
%% Generate sensitivity plot (centre of fascicle 1) 

  clf, hold on, m = mdl_(f_id);
  sel = m.object_id{strncmpi(m.object_name,'PDMS',4)}; 
  idx = m.elems(sel,:);
  idx = idx(sum(reshape(m.nodes(idx,2), [], 4) >= -10*eps,2) > 2,:);

  y_gscale = [1 .25 .35];
  if min(EM.(fascicleN).pot(:)) < 0.5  
    y_gscale(3) = max(EM.(fascicleN).pot(:)) + 0.025 - min(m.nodes(idx,3))*y_gscale(2);
  end

  patch('Faces',idx, 'Vertices',m.nodes(:,[1 3]) .* y_gscale(1:2) + [0 y_gscale(3)], ...
        'EdgeColor',[.7 .7 .7], 'FaceColor','w','EdgeAlpha',0.2)

  outline = m.nodes(unique(idx),[1 3]);
  idx = convhull(outline); 
  
  plot(outline(idx,1),outline(idx,2)*y_gscale(2) + y_gscale(3), '-', ...
      'Color', [.5 .5 .5],'LineWidth',1.2)
    
  if isfield(EM.info,'FascicleTrajectory')    
      alist = dir(tools.file('in~\axons\*.mat'));    
      error get_fascicle_Trajectory_anatomy
      % load(tools.INPUT_file(alist,file),'F');
      
      xy = mean(F.fascicles); 
      xyz = tools.from_trajectory(EM,F,xy);
      xyz = permute(xyz(:,[3 2 1],:),[3 2 1]);
      z = cumsum(sqrt(sum(diff(xyz).^2,2))); % length of axon
      z(end+1) = [2 -1]*z([end end-1]); % extend by 1 (because diff)
      z = z - z(find(xyz(:,3) == min(xyz(:,3)),1));
  else   
      xy = [mean(z_(f_id))   mean(y_(f_id))  0];  
      z = linspace(-8,8,401);
      xyz = xy + z' * [0 0 1];
  end
  
  for ee = 1:nE

      pot = elec_sensor{ee}(xyz);

      plot(z,pot,'Linewidth',1.2,'Color',C(ee,:)), hold on

      e_n = m.electrode(ee).nodes;
      e_xy = [min(m.nodes(e_n,1)) min(m.nodes(e_n,3))*y_gscale(2) + y_gscale(3) ... 
              range(m.nodes(e_n,1)) range(m.nodes(e_n,3))*y_gscale(2)]; 

      e_circle = 0; 
      rectangle('Position',e_xy,'FaceColor',W(ee,1),'EdgeColor',C(ee,:), ... 
                'Curvature',[1 1]*e_circle,'LineWidth',1.2)
  end

  ok = ~isnan(pot); 
  plot(z(ok),xyz(ok,1)*y_gscale(2) + y_gscale(3),'--','Color',G(3),'LineWidth',1.2)

  set(gca,'DataAspectRatio',y_gscale./min(y_gscale))
  axis tight, tools.tidyPlot
  xlabel('Distance in space, mm')
  ylabel('Sensitivity, �V/�A')
  [~,lbl] = fileparts(EM.filename);
  title(strrep(lbl,'_','\_'))

  clear xyz xy z y_gscale sel idx pot ee e_n e_xy ans outline
return

%% Cross-section of fascicles
function choose_sensitivity_figure(EM, opts) 

if nargin < 2, opts = struct; end
if ~isfield(opts,'fascicle'),  opts.fascicle = 1; end
if ~isfield(opts,'electrode'), opts.electrode = []; end
if ~isfield(opts,'xy'),        opts.xy = []; end
if ~isfield(opts,'do_hwidth'), opts.do_hwidth = true; end
if ~isfield(opts,'do_parametric'), opts.do_parametric = false; end
if ~isfield(opts,'do_image'),  opts.do_image = false; end
if ~isfield(opts,'force_mf'),  opts.force_mf = false; end
if ~isfield(opts,'do_PCA'),    opts.do_PCA = 0; end
if ~isfield(opts,'pk_quantile'), opts.pk_quantile = 1; end

if nargin < 1 && evalin('caller','exist(''EM'',''var'')')
    EM = evalin('caller','EM');    
end

nE = sum(~cellfun(@isempty,{EM.model(1).electrode.name}));
if isempty(opts.electrode), opts.electrode = 1:nE; end % default: all elec

if numel(opts.fascicle) > 1 || opts.force_mf, % Generate multiple-fascicle figure
    make_panels_multiFascicle(EM,opts)
    return
elseif numel(EM.model) > 1, EM.model = EM.model(opts.fascicle);
end

% Generate interactive figure(s)
if numel(opts.electrode) == 1, make_panels_debugMeshInterp(EM,opts);
else                           make_panels_singleFascicle(EM,opts);
end
return

%% Interactive (clickable) figure for single fascicle anatomy
function make_panels_singleFascicle(EM, opts)

e_id = opts.electrode;
f_id = opts.fascicle;

nE = sum(~cellfun(@isempty,{EM.model(1).electrode.name}));

% Standard color-tools for figures
C = lines(nE); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

% Create utility local functions

% for f = fieldnames(EM.utils)' 
%     EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
%     eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
% end

x_ = @(i) EM.model.nodes(i,1);
y_ = @(i) EM.model.nodes(i,2);
z_ = @(i) EM.model.nodes(i,3);
fac_ = @(i) EM.(sprintf('Fascicle%d',i));

if evalin('caller','exist(''sensor'',''var'')') % might be preloaded.
  sensor = evalin('caller','elec_sensor');      % get from caller.
else                                            % Generate if not. 
  sensor(e_id) = {[]}; % allocate cell array
  if opts.debug_thinlayer_idx % remove P_Fasc from interpolation
       idx = (strcmp(EM.model.object_name,sprintf('Fascicle%d',f_id)));
       idx = unique(EM.model.elems(EM.model.object_id{idx},:));
  else idx = fac_(f_id).idx; 
  end 
  fok = ismember(fac_(f_id).idx,idx(:)); 
  idx = fac_(f_id).idx(fok); 
  
  sensor{e_id(1)} = scatteredInterpolant(z_(idx), y_(idx), x_(idx), ...
                                      fac_(f_id).pot(e_id(1), fok)', ...
                                                'natural','none');
  for ee = e_id 
     sensor{ee} = sensor{e_id(1)}; 
     sensor{ee}.Values = fac_(f_id).pot(ee, fok)';
   end
   assignin('caller','elec_sensor',sensor)
end

nE = min(nE, numel(e_id));
%%

xy = []; 
a_list = dir(tools.file('axons~/ax*.mat')); 

nP = size(opts.xy,1);
if nP == 1, nP = opts.xy; elseif nP > 1, xy = opts.xy; end

if isempty(xy) && numel(a_list) > 0
  %% Load from file
  ax_file = tools.file('get','axons~/ax*.mat','newest'); 
  ax = load(ax_file);
  if isfield(ax,'F')
    if isfield(ax.F,'fascicles')
         xy_fac = ax.F.fascicles(:,:,f_id); % n x 2
    else xy_fac = ax.F.outline(:,:,f_id); % n x 
    end
    
    if max(abs(xy_fac(:))) > max(abs(sensor{1}.Points(:)))
      xy_fac = xy_fac / 1e3; % um to mm
    end
    xy = [ax.pop.fibre_xy(ax.pop.fibre_fascicle == f_id,:);...
          ax.pop.unmyelinated_xy(ax.pop.unmyelinated_fascicle == f_id,:)]/1e3;

     [~,sel] = sort(rand(size(xy(:,1))));
     xy = xy(sel(1:min(400,end)),:);
     
     ok = in_loop(xy,xy_fac); 
     if mean(ok) < 0.9, 
       
       warning('axons found in %s did not match outline specified in %s', ...
                strrep(ax_file,tools.file,'~'),ax.F.filename)
       xy = []; 
     else
       fprintf('using axons found in %s\n',strrep(ax_file,tools.file,'~'))
     end
  end
end
if isempty(xy)
  
    %% fake up the fascicle outlines, axon xy

  xy_fac = [z_(fac_(f_id).idx) y_(fac_(f_id).idx)];
  sel = convhull(xy_fac(:,1),xy_fac(:,2));
  if numel(sel) > 100, sel = sel(round(linspace(1,end,101))); end
  xy_fac = xy_fac(sel([1:end 1]),:);
  
  
  if nP == 0
    
    xy0 = [median(xy_fac) range(xy_fac)];

    [gx,gy] = meshgrid(xy0(1)-xy0(3)*linspace(-1,1,opts.xy_res), ...
                       xy0(2)-xy0(4)*linspace(-1,1,opts.xy_res));

    if isfield(EM.info,'FascicleTrajectory')

      xyz = tools.from_trajectory(EM,F,[gx(:) gy(:)]);
      ok = in_loop([gx(:) gy(:)],F.outline([1:end 1],:,ff));
      ok = ok & ~isnan(sensor{ff}(xyz(:,[3 2 1],round(end/2))));
    else 
      ok = in_loop([gx(:) gy(:)],xy_fac);
      if ~any(ok), ok = in_loop([gx(:) gy(:)],flipud(xy_fac));end
    end
    xy = [gx(ok) gy(ok)]; 
  else
    xy = rand(nP,2) .* range(xy_fac) + min(xy_fac);
    ok = in_loop(xy,xy_fac); 
    xy = xy(ok,:); 
  end


end

%%

nP = size(xy,1); 
z = linspace(-6,6,201);
% z = [0:0.05:6]; z = [-fliplr(z) z(2:end)];
dz = mean(diff(z)); 

halfwidth = zeros(nP,nE,1);
peak_val  = zeros(nP,nE,1);
edge_val  = zeros(nP,nE,2);

%% Right side: per-electode view

h_elec = gobjects(0); 

printInfo(); 

for ee = 1:nE
    %%
    h_elec(ee) = subplot(nE,2,2*ee); hold on
    
    if opts.do_PCA > 0, potential = zeros(numel(z),size(xy,1)); end
    
    for ii = 1:size(xy,1)

        printInfo('elec%d : %0.2f%% ... ', e_id(ee),100*ii/size(xy,1))
        xyz = [xy(ii,:) 0] + z' * [0 0 1];
        pot = sensor{e_id(ee)}(xyz);
        
        if all(isnan(pot)), peak_val(ii,ee) = nan; continue, end

        if opts.do_PCA == 0 
          plot(z,pot,'-','Color',[W(ee,8*(xy(ii,2) - min(xy(:,2)))) 0.5],'UserData',ii)
        else
          fix = isnan(pot); 
          if any(fix) && mean(fix) < 0.5
            pot(fix) = interp1(z(~fix),pot(~fix),z(fix),'linear','extrap');
          end
          potential(:,ii) = pot; 
        end
        
        pot(isnan(pot)) = [];
        
        pv = max(pot);
        ev = pot([1 end]);
        
        pot = abs(pot) - min(abs(pot)); 
        
        peak_val(ii,ee)   = pv;
        edge_val(ii,ee,:) = ev;
        halfwidth(ii,ee) = dz * range(find(pot > max(pot)/2));         
    end
    
    if opts.do_PCA > 0
    
      
      [coeff,profile,latent] = pca(potential);
      scale = max(abs(coeff)); 
      latent = latent / sum(latent);
      
      profile = profile .* (scale);
      coeff   = coeff ./ (scale);
      
      for kk = 1:opts.do_PCA
         plot(z,profile(:,kk),'-','Color',C(kk,:),'LineWidth',1.1,'UserData',latent(kk))
      end
      coeff = coeff(:,1:opts.do_PCA);
      
    else
      coeff = [peak_val(:,ee) halfwidth(:,ee)];
    end
    
    axis tight, tools.tidyPlot
    xlim([-5.1 5]), set(gca,'XTick',-5:5,'UserData',ee);    
    if ee < nE, set(gca,'XTickLabel',''),
    else xlabel('Distance, mm')
    end
    
    set(gca,'Color','w','UserData',coeff, ...
            'ButtonDownFcn',@(a,b) update_anatomy_plots(a,b,xy_fac,xy))
          
end

printInfo(); fprintf('Done!\n')

set(h_elec,'YLim',[min([h_elec.YLim]) max([h_elec.YLim])])

[~,t] = fileparts(EM.filename);
t = regexprep(t,'([_\\])','\\$1');

h = suptitle(t); % ,'fontsize',9)
h.FontSize = 12; h.Color = [.5 .5 .5];

update_anatomy_plots(gca,[],xy_fac,xy)

    

function update_anatomy_plots(hobj,~,outline,xy)


coeff = hobj.UserData;
bad = any(isnan(coeff),2);

xy(bad,:) = []; 
coeff(bad,:) = [];

nY = size(coeff,2);

if nY == numel(hobj.Children)
     labels = arrayfun(@(n) sprintf('PCA-%d (au)',n), 1:nY,'unif',0);   
elseif nY == 2, labels = {'Peak, �V/�A','Half-width, mm'};
else labels = arrayfun(@(n) sprintf('y%d (?)',n), 1:nY,'unif',0);
end

%% Set up left hand side axes

for yy = 1:nY

  h = subplot(nY,2,2*yy-1); cla reset
  plot(outline(:,1),outline(:,2),'-','Color',[0 0 0 0.3])
  axis equal, hold on

  ylim([0 max(ylim)]), tools.tidyPlot
  scatter(xy(:,1),xy(:,2),120,coeff(:,yy),'.')
  ch = colorbar; ylabel(ch,labels{yy}); 

  if yy == 1, set(gca,'UserData','peak'), end

% rectangle('Position', [0 min(y_(EM.model.electrode(1).nodes)) max(xlim) ...
%                        range(y_(EM.model.electrode(1).nodes))], ...
%           'FaceColor', G(4), 'EdgeColor','none')
end

colormap magma

%% Do something different for multple fascicles
function make_panels_multiFascicle(EM, opts)


% f_id, e_id, xy

nF = max(opts.fascicle); 
nE = numel(opts.electrode);
nP = size(opts.xy,1);
if nP == 1, nP = opts.xy; end

FLAG_multimesh = (numel(EM.model) > 1); 
FLAG_stimulus = isfield(EM,'v_extracellular');

% % Create utility local functions
% for f = fieldnames(EM.utils)' 
%     EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
%     eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
% end

fac_ = @(i)EM.(sprintf('Fascicle%d',i)).idx;

if FLAG_multimesh
  x_ = @(f,i) EM.model(f).nodes(i,1);
  y_ = @(f,i) EM.model(f).nodes(i,2);
  z_ = @(f,i) EM.model(f).nodes(i,3); 
else
  x_ = @(f,i) EM.model(1).nodes(i,1);
  y_ = @(f,i) EM.model(1).nodes(i,2);
  z_ = @(f,i) EM.model(1).nodes(i,3);
end

% Standard color-tools for figures
% C = lines(max(f_id)); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

%%

clf
z = linspace(-5,5,101);
dz = mean(diff(z)); 
xy = cell(nF,1); 
outline = cell(nF,1); 

sensor = cell(nF,1); 

for ee = 1:nE
      
  fwhm_val = cell(nF,1);
  peak_val  = cell(nF,1);
  edge_val  = cell(nF,1);
    
  for f0 = 1:numel(opts.fascicle)
    
    ff = opts.fascicle(f0); 
      
    e_pot = EM.(sprintf('Fascicle%d',ff)).pot(opts.electrode(ee),:)';
    if isempty(sensor{ff})
         sensor{ff} = scatteredInterpolant(z_(ff,fac_(ff)), y_(ff,fac_(ff)), ... 
                                           x_(ff,fac_(ff)), ...
                                           e_pot,'natural','none'); 
    else sensor{ff}.Values = e_pot;
    end
    
    is_3D_trajectory = isfield(EM,'info') && ...
                       isfield(EM.info,'FascicleTrajectory');
    
    if isempty(xy{ff})
      if is_3D_trajectory
      
        alist = dir(tools.file('in~\axons\*.mat'));
        error get_FascicleTrajectory
        % load(tools.INPUT_file(alist,EM.filename),'F');
        
        xyf = F.outline(:,:,ff);
        outline{ff} = F.outline(:,:,ff);
      else      
        xyf = [z_(ff,fac_(ff)) y_(ff,fac_(ff))];
        idx = convhull(xyf(:,1),xyf(:,2));
        outline{ff} = xyf(idx([1:end 1]),:);
      end

      if nP == 0 || opts.do_image

        xy0 = [median(xyf) min(xyf) max(xyf)];
        xy0(3:4) = max(abs(xy0(1:2)-xy0([3 4; 5 6])),[],1);

        [gx,gy] = meshgrid(xy0(1)-xy0(3)*linspace(-1,1,opts.xy_res), ...
                           xy0(2)-xy0(4)*linspace(-1,1,opts.xy_res));
                         
                         
        if is_3D_trajectory
          
          xyz = tools.from_trajectory(EM,F,[gx(:) gy(:)]);
          ok = in_loop([gx(:) gy(:)],F.outline([1:end 1],:,ff));
          ok = ok & ~isnan(sensor{ff}(xyz(:,[3 2 1],round(end/2))));
        elseif opts.do_image, ok = true(size(gx(:))); 
        else ok = ~isnan(sensor{ff}([gx(:) gy(:)*[1 0]]));
        end
                         
        xy{ff} = [gx(ok) gy(ok)]; 
      else
        xy{ff} = rand(ceil(nP*4/pi),2) .* range(xyf) + min(xyf);
      end
      
    end

    peak_val{ff} = zeros(size(xy(:,1)));
    fwhm_val{ff} = zeros(size(xy(:,1)));
    edge_val{ff} = zeros(size(xy));    
    
    for ii = 1:size(xy{ff},1)
        
        if is_3D_trajectory
          
          
          xyz = tools.from_trajectory(EM,F,xy{ff}(ii,:));
          xyz = permute(xyz(:,[3 2 1],:),[3 2 1]);
          z = cumsum(sqrt(sum(diff(xyz).^2,2))); % length of axon
          z(end+1) = [2 -1]*z([end end-1]); % extend by 1 (because diff)
          z = z - z(find(xyz(:,3) == min(xyz(:,3)),1));
          dz = mean(diff(z));
        else
          xyz = [xy{ff}(ii,:) 0] + z' * [0 0 1];
        end
        
        
        pot = sensor{ff}(xyz);
        
        if mean(isnan(pot)) > 0.9, peak_val{ff}(ii,1) = nan; continue, end

        % plot(z,pot,'-'),'Color',[W(ff,8*(xy{ff}(ii,2) - min(xy{ff}(:,2)))) 0.5])

        pot(isnan(pot)) = [];        
        [pv,x0] = max(pot);
        x0 = z(x0); 
        
        if opts.pk_quantile < 1
          pv = quantile(pot,opts.pk_quantile);
        end

        ev = pot([1 end]);

        pot = abs(pot) - min(abs(pot)); 

        peak_val{ff}(ii,1)   = pv;
        edge_val{ff}(ii,1:2) = ev;
        fwhm_val{ff}(ii,1) = dz * range(find(pot > max(pot)/2));   
        
        if opts.do_parametric % I can't recommend this, fit quality isn't great
          
          x_skew = @(x,x0,w) (x-x0).*(1+w*(x>x0)); 
          Cauchy_fit = @(x,p) p(1).*p(2) ./ ( x_skew(x,p(5),p(4)).^2 + ...
                                              + p(2).^2) / pi + p(3); 
          lsq = @(p) mean((pot' - Cauchy_fit(z,[p x0])).^2);
          
          b0 = quantile(pot,0.1); 
          
          pInit = [2*pv fwhm_val(ii,ff) b0 0]; 
          
          LB = [0 0 0 -1];
          UB = [inf range(z) max(pot) 1];
          
          fopts = optimset('Display','off');          
          pFit = fmincon(lsq,pInit,[],[],[],[],LB,UB,[],fopts);
          
          peak_val{ff}(ii,1) = pFit(1);
          fwhm_val{ff}(ii,1) = pFit(2);          
          edge_val{ff}(ii,1:2) = pFit(3:4);           
          %%
          % cla, plot(z,pot,'.'), hold on          
          % plot(z,Cauchy_fit(z,[pFit x0]))
          % plot(z,Cauchy_fit(z,[pInit x0]),'--')
        end
    end
    
    if opts.do_hwidth, subplot(nE,2,2*ee-1),
    else               subplot(ceil(nE/2),2,ee),
    end,               hold on, axis equal
    
    if opts.do_image, peak_val{ff} = reshape(peak_val{ff},size(gx,1),[]);
      imagesc(gx(1,:),gy(:,1),peak_val{ff})
      plot(outline{ff}(:,1),outline{ff}(:,2),'-','Color',[0 0 0 0.6])
      colormap(gca,[.9 .9 .9; magma(256)])
    else
      plot(outline{ff}(:,1),outline{ff}(:,2),'-','Color',[0 0 0 0.3])
      scatter(xy{ff}(:,1),xy{ff}(:,2),[],peak_val{ff}(:),'.')
      colormap(gca,tools.magma)
    end
    
    fwhm_val{ff}(isnan(peak_val{ff})) = NaN;
    
    if ~opts.do_hwidth, continue, end
    
    subplot(nE,2,2*ee), hold on, axis equal
    
    if opts.do_image, fwhm_val{ff} = reshape(fwhm_val{ff},size(gx,1),[]);
      imagesc(gx(1,:),gy(:,1),fwhm_val{ff})
      plot(outline{ff}(:,1),outline{ff}(:,2),'-','Color',[0 0 0 0.6])
      colormap(gca,[.9 .9 .9; tools.magma(256)])
    else
      plot(outline{ff}(:,1),outline{ff}(:,2),'-','Color',[0 0 0 0.3])
      scatter(xy{ff}(:,1),xy{ff}(:,2),[],fwhm_val{ff}(:),'.')
      colormap(gca,tools.magma)
    end
  end
  
  %% Format axes
  
  if opts.do_hwidth, subplot(nE,2,2*ee-1)
  else               subplot(ceil(nE/2),2,ee),
  end
  
  axis tight, ylim([0 max(ylim)]), tools.tidyPlot
  caxis(quantile(cat(1,peak_val{:}),[.01 .99]))  
  
  if FLAG_stimulus, % generate correct stim label
    stim = EM.model(1).stimulation(opts.electrode(ee));
    lbl = sprintf('\\bfS%s,R%s,',sprintf('%d+',find(stim.stim_pattern == 1)), ...
                             sprintf('%d+',find(stim.stim_pattern == -1)));
       ylabel(strrep(lbl,'+,',''))
  else ylabel(sprintf('\\bfE%d',opts.electrode(ee)))  
  end
  
  if ee == 1, title('Peak (�V/�A)'), end

  if ~opts.do_hwidth, pause(0.02), continue, end
  
  ch = colorbar; ch.Position(1) = 0.424;
  
  subplot(nE,2,2*ee), ylim([0 max(ylim)]), tools.tidyPlot
  vals = quantile(cat(1,fwhm_val{:}),[.01 .99]);
  if numel(unique(vals)) == 1 || any(isnan(vals))
    vals = [nanmin(cat(1,fwhm_val{:}))-1e-3 ...
            nanmax(cat(1,fwhm_val{:}))+1e-3];
  end
  caxis(vals)
  ch = colorbar; ch.Position(1) = 0.904;
  if ee == 1, title('Half-width (mm)'), end
  
  pause(0.02);
  
end

return

%% Utility Functions
function EM = convert_stim2fascicles(EM)
  
  FLAG_multimesh = (numel(EM.model) > 1); 
  named = evalin('caller','named');
  

  if FLAG_multimesh, EM.model = [EM.model{:}];
       nF = numel(EM.model); 
  else nF = sum(contains(EM.model.object_name,'Fasc')) - ...
            sum(contains(EM.model.object_name,'P_Fasc'));
  end
  
  fasc_ = @(n) sprintf('Fascicle%d',n);
  
  % EM.v_extracellular is a list, (nodes x electrodes) of the sequential
  % bipolar stimulus-induced extracellular potential, for ALL nodes 
  % (not just fascicle nodes). 

  % a stimulus file is broken up by fascicle, which each fascicle
  % containing a NODE index and a list of per-node values 
  
  % EM.model.object_id is a list of ELEMENT ids which make up each object
  % in the EIDORS model. Gather each NODE in that set of elements, then use
  % that to emulate the EM.FascicleN.pot/idx 
    
  for ff = 1:nF
    
    if FLAG_multimesh, m = ff; else m = 1; end
    
    sel = strcmpi(EM.model(m).object_name,fasc_(ff));
    if sum(sel) == 0 && FLAG_multimesh, 
      sel = strcmpi(EM.model(m).object_name,fasc_(1));      
    end
    if sum(sel) ~= 1, error('Object %s not found in EM.model.object_name',fasc_(ff)); end
    idx = unique(EM.model(m).elems(EM.model(m).object_id{sel},:));
    
    EM.(fasc_(ff)).idx = idx'; 
    if FLAG_multimesh
         EM.(fasc_(ff)).pot = EM.v_extracellular{ff}(idx,:)';
    else EM.(fasc_(ff)).pot = EM.v_extracellular(idx,:)';
    end
  end
  
  if any(named('-AF'))
    mk_shortcuts(EM,FLAG_multimesh)
    for ff = 1:nF % ACTIVATING FUNCTION (if requested) 
      dx = 1e-3; 
      Ve = scatteredInterpolant(z_(ff), y_(ff), x_(ff), ...
                                EM.(fasc_(ff)).pot(1,:)','linear','linear'); 
      for ee = 1:nE

        Ve.Values = EM.(fasc_(ff)).pot(ee,:)';
        xyz = [z_(ff) y_(ff) x_(ff)-dx];  vx0 = Ve(xyz); % have to 2-step
        xyz = [z_(ff) y_(ff) x_(ff)+dx];  vx1 = Ve(xyz);

        EM.(fasc_(ff)).pot(ee,:) = (vx0+vx1-2*Ve.Values)/(dx^2);
      end
    end
  end
  
return
function EM = convert_mono2bipolar(EM)

  FLAG_multimesh = (numel(EM.model) > 1);  

  if FLAG_multimesh, EM.model = [EM.model{:}];
       nF = numel(EM.model); 
  else nF = sum(contains(EM.model.object_name,'Fasc')) - ...
            sum(contains(EM.model.object_name,'P_Fasc'));
  end
  
  nE = size(EM.Fascicle1.pot,1);
  fasc_ = @(n) sprintf('Fascicle%d',n);
  ret_ = @(e) mod(e,nE)+1; % sequentual bipolar pattern

  for ff = 1:nF, mono = EM.(fasc_(ff)).pot;
    for ee = 1:nE
      EM.(fasc_(ff)).pot(ee,:) = mono(ee,:) - mono(ret_(ee),:); 
    end
  end
  for ee = 1:nE % set "stim" structure 
    
    this = struct;
    this.stim_pattern = zeros(nE,1);
    this.stim_pattern(ee) = 1;
    this.stim_pattern(ret_(ee)) = -1;   
    if isempty(EM.model.stimulation)
         EM.model.stimulation = this; 
    else EM.model.stimulation(ee) = this;
    end
    
  end
  
  EM.v_extracellular = false; 
  
return


function mk_shortcuts(EM,is_multimesh)

if evalin('caller','exist(''x_'',''var'')'), return, end
if nargin == 1, is_multimesh = (numel(EM.model) > 1); end

if is_multimesh
  
  mdl_ = @(n) EM.model(n);
  fac_ =  @(n) EM.(sprintf('Fascicle%d',n)).idx;
  x_ = @(n) EM.model(n).nodes(fac_(n),1); 
  y_ = @(n) EM.model(n).nodes(fac_(n),2); 
  z_ = @(n) EM.model(n).nodes(fac_(n),3);
else  
  mdl_ = @(n) EM.model(1);
  fac_ =  @(n) EM.(sprintf('Fascicle%d',n)).idx; 
  x_ = @(n) EM.model.nodes(fac_(n),1); 
  y_ = @(n) EM.model.nodes(fac_(n),2); 
  z_ = @(n) EM.model.nodes(fac_(n),3); 
    
  % for f = fieldnames(EM.utils)' % Create utility local functions
  %   EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
  %   eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
  % end
end

nE = size(EM.Fascicle1.pot,1);
nF = sum(contains(fieldnames(EM),'Fascicle'));

% Standard color-tools for figures
C = lines(max(7,nE)); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

for k = {'mdl_','fac_','x_','y_','z_','nE','nF','C','W','G'}
  assignin('caller',k{1},eval(k{1}));
end
function is = in_loop(xy,loop) % are xy contained in the loop? 

is = false(size(xy(:,1)));

for ii = 1:size(xy,1)
    
    dx = sign(loop(:,1) - xy(ii,1));
    cx = find(dx ~= circshift(dx,[1 0])); 
    if isempty(cx), continue, end
    cn = 0*cx; 
        
    % clf, hold on, C = lines(7);
    % plot(path(:,1),path(:,2),'color',[0 0 0 0.3])
    % plot(xy(ii,1),xy(ii,2),'s','LineWidth',1.2,'Color',C(2,:))
    % plot(path(cx,1),path(cx,2),'.k','MarkerSize',8)    
    % plot([1 1]*xy(ii,1),ylim,'--','LineWidth',1.2,'Color',C(2,:))
    
    for pp = reshape(cx,1,[])        
        % plot(path(pp-[0 1],1),path(pp-[0 1],2),'-k','LineWidth',1.2)
        
        seg = loop(mod(pp-[1 2],size(loop,1))+1,:);
        
        if all(seg(:,2) >= xy(ii,2)), cn(cx==pp)=1; continue
        elseif all(seg(:,2) < xy(ii,2)),            continue
        end
        
        u = interp1(seg(:,1),seg(:,2),xy(ii,1));
        cn(cx==pp) = (u > xy(ii,2));
    end
    
    is(ii) = mod(sum(cn),2)==1;
end




%% Interactive DEBUG mesh viewer
function make_panels_debugMeshInterp(EM,opts)

% related to make_singleFascicle_panels
e_id = opts.electrode;
f_id = opts.fascicle;

nE = sum(~cellfun(@isempty,{EM.model(1).electrode.name}));
nE = min(nE, numel(e_id));

% Standard color-tools for figures
C = lines(7); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

% Create utility local functions
x_ = @(i) EM.model.nodes(i,1);
y_ = @(i) EM.model.nodes(i,2);
z_ = @(i) EM.model.nodes(i,3);
fac_ = @(i) EM.(sprintf('Fascicle%d',i));

if evalin('caller','exist(''sensor'',''var'')') % might be preloaded.
  sensor = evalin('caller','elec_sensor');      % get from caller.
else                                            % Generate if not. 
  sensor(e_id) = {[]}; % allocate cell array
  if opts.debug_thinlayer_idx % remove P_Fasc from interpolation
       idx = (strcmp(EM.model.object_name,sprintf('Fascicle%d',f_id)));
       idx = EM.model.elems(EM.model.object_id{idx});
  else idx = fac_(f_id).idx; 
  end 
  fok = ismember(fac_(f_id).idx,idx(:)); 
  idx = fac_(f_id).idx(fok); 
  
  sensor{e_id(1)} = scatteredInterpolant(z_(idx), y_(idx), x_(idx), ...
                                      fac_(f_id).pot(e_id(1), fok)', ...
                                                'natural','nearest');
  for ee = e_id 
     sensor{ee} = sensor{e_id(1)}; 
     sensor{ee}.Values = fac_(f_id).pot(ee, fok)';
   end
   assignin('caller','elec_sensor',sensor)
end

%% Get fascicle and perineurium outline

nerve_xy = cell(1,3); 
for layer = 1:3 % all FASC, bndry FASC, bndry P_FASC
  
  switch layer
    case 1
      idx = (strcmp(EM.model.object_name,sprintf('Fascicle%d',f_id)));
      fascicle_indices = EM.model.object_id{idx}; 
      idx = unique(EM.model.elems(fascicle_indices,:));
      
    case 2
      idx = unique(EM.model.boundary(ismember(EM.model.boundary(:),idx(:))));
    case 3
      idx = (strcmp(EM.model.object_name,sprintf('P_Fascicle%d',f_id)));
      if ~any(idx), break, end
      idx = unique(EM.model.elems(EM.model.object_id{idx},:));

      % idx = unique(EM.model.boundary(ismember(EM.model.boundary(:),idx(:))));
  end
 
  nerve_xy{layer} = [z_(idx) y_(idx)];
  sel = convhull(nerve_xy{layer}(:,1),nerve_xy{layer}(:,2));
  if numel(sel) > 100, sel = sel(round(linspace(1,end,101))); end
  nerve_xy{layer} = nerve_xy{layer}(sel([1:end 1]),:);
end


xy = median(nerve_xy{2}); 
z = (0:0.02:6); z = [-fliplr(z) z(2:end)];

%% Panel 1: fascicle cross-section (clickable)
clf
subplot(2,2,1), cla, hold on

if ~isempty(nerve_xy{3})
  plot(nerve_xy{3}(:,1),nerve_xy{3}(:,2),'-','Color',G(1),'LineWidth',1)
end

plot(nerve_xy{2}(:,1),nerve_xy{2}(:,2),'-','Color',G(6),'LineWidth',1.3)
plot(nerve_xy{1}(:,1),nerve_xy{1}(:,2),'-','Color',[0 0 0 0.3])

plot(xy(1),xy(2),'o','Color',G(6),'LineWidth',1.3)
plot(xy(1),xy(2),'+','Color',C(2,:),'LineWidth',1,'MarkerSize',12)

axis image xy, tools.tidyPlot
set(gca,'Color','w','ButtonDownFcn',@update_mesh_view);
set(get(gca,'Children'),'hittest','off')
% set(gca,'UserData',fascicle_indices); 
title('Click to set XY')

%%
subplot(2,2,2), cla, hold on

xyz = [xy 0] + z' * [0 0 1];
pot = sensor{e_id(1)}(xyz);

plot(z,pot,'-','Color',G(6),'LineWidth',1.3,'UserData','median')
plot(z,pot,'-','Color',[C(2,:) 0.7],'LineWidth',1,'UserData','active')

tools.tidyPlot
set(gca,'UserData',sensor{e_id(1)})
set(gca,'Color','w','ButtonDownFcn',@update_mesh_axes);
set(get(gca,'Children'),'hittest','off')
title('Click to set ROI')

init.Button = 0;
init.IntersectionPoint = [xy 0]; 
%%
subplot(2,2,[3 4]), cla, hold on

% idx = unique(pointLocation(h(3).UserData,xyz)); % (:,[3 2 1])));
tets = EM.model.elems(fascicle_indices,[1 2 3 2 4 1 1 3 4]);

px = reshape(EM.model.nodes(tets',1),3,[]); px(4,:) = NaN;
py = reshape(EM.model.nodes(tets',2),3,[]); py(4,:) = NaN;
pz = reshape(EM.model.nodes(tets',3),3,[]); pz(4,:) = NaN;

plot3(px(:),pz(:),py(:),'Color',[0 0 0 0.05])
plot3(xyz(:,3),xyz(:,2),xyz(:,1),'Color',[0 0 0 0.5])

scatter3(xyz(:,3),xyz(:,2),xyz(:,1),10,pot,'o','filled')
scatter3(xyz(1,3),xyz(1,2),xyz(1,1),8,pot(1),'s')
axis image, tools.tidyPlot, grid on, colormap(tools.magma)
% set(get(gca,'Children'),'Clipping','off')
set(gca,'UserData',triangulation(double(EM.model.elems),EM.model.nodes));
update_mesh_view(gca,init);

xlabel('x'), ylabel('y'), zlabel('z')

%% Panel 2: trace along mesh
function update_mesh_view(h,e)

h = flipud(h.Parent.Children);
%%

z = h(2).Children(1).XData;
xy  = e.IntersectionPoint(1:2);
xyz = [xy 0] + z' * [0 0 1];
pot = h(2).UserData(xyz);
h(2).Children(1).YData = pot;
h(1).Children(1).XData = xy(1);
h(1).Children(1).YData = xy(2);
%%

% h(1) is the node_index list for fascicleN
% h(2) is the scatteredInterpolant for fascicleN
% h(3) is a triangulation from EM.model

idx = unique(pointLocation(h(3).UserData,fliplr(xyz)));
idx(isnan(idx)) = []; 

% idx = unique(pointLocation(h(3).UserData,xyz)); % (:,[3 2 1])));
tets = h(3).UserData.ConnectivityList(idx,[1 2 3 2 4 1 1 3 4]);

px = reshape(h(3).UserData.Points(tets',1),3,[]); px(4,:) = NaN;
py = reshape(h(3).UserData.Points(tets',2),3,[]); py(4,:) = NaN;
pz = reshape(h(3).UserData.Points(tets',3),3,[]); pz(4,:) = NaN;

h(3).Children(3).XData = px(:);
h(3).Children(3).YData = pz(:);
h(3).Children(3).ZData = py(:);

h(3).Children(2).XData = xyz(:,3);
h(3).Children(2).YData = xyz(:,1);
h(3).Children(2).ZData = xyz(:,2);
h(3).Children(2).CData = pot;

xyz = h(3).UserData.Points(unique(tets(:)),:);
pot = h(2).UserData(fliplr(xyz));

h(3).Children(1).XData = xyz(:,1);
h(3).Children(1).YData = xyz(:,3);
h(3).Children(1).ZData = xyz(:,2);
h(3).Children(1).CData = pot;

% not_fasc = unique(tets(ismember(tets(:),h(1).UserData)));
% both tets and h(1).UserData /should/ be node indices, so IDK why I
% couldn't use ismember[] to pull out the perineurium node indices

%%
return
function update_mesh_axes(h,~)


h = flipud(h.Parent.Children);

xy = [h(1).Children(1).XData h(1).Children(1).YData];


expand = [1 0; 0 1] + [1 -1; -1 1]*1;

caxis(h(3),ylim(h(2)))
xlim(h(3),xlim(h(2)))
ylim(h(3),xlim(h(1))*expand)
zlim(h(3),ylim(h(1))*expand)

return

