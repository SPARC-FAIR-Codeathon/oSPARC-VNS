
function stats = sequential_F_test(X, Y, varargin)
% Statistical approach used in Eiber 2021a 

named = @(v) strncmpi(v,varargin,length(v));
get_ = @(v) varargin{find(named(v))+1};

if nargin < 2, error('Not enough input arguments'), end

do_PDF = any(named('-pdf')); 
do_LOGSCALE = any(named('-log')); 

if any(named('-boot')) && isstruct(X)
  stats = bootstrap_parameters(Y, varargin{1}, X);
  return
end

if do_PDF, ps_folder = 'make-pdf\';   
  if any(named('-pdf-chunk')), tt = get_('-pdf-chunk');
   if ~exist([tempdir ps_folder],'dir'), mkdir([tempdir ps_folder]); end
  else tt = 0;
   if exist([tempdir ps_folder],'dir'), rmdir([tempdir ps_folder],'s'); 
   end,                                 mkdir([tempdir ps_folder]);
  end
end

[equations, labels, n_terms] = setup_statistical_models; 

if any(named('-optimset')), f_opt = get_('-optimset'); 
  if iscell(f_opt), f_opt = optimset(f_opt{:}); end
else f_opt = optimset('maxFunEvals',1e6,'maxiter',1e6);
end

%% For each axon class run the above battery of equations and produce results

Y_label = 'Y';
if any(named('-title')), Y_label = get_('-title'); end
fprintf('\n running F TEST battery: %s \n', Y_label)

if do_LOGSCALE, Y = log10(Y);
  if any(Y_label == '('), Y_label = strrep(Y_label,'(','(log_{10} ');
  else Y_label = ['log_{10} ' Y_label]; 
  end
end

ok = ~isnan(Y); Y = Y(ok); X = X(ok,:);

%% Fit models and compute residuals 

isa_cat_var = cellfun(@(x) all(ismember(x,min(x):max(x))), num2cell(X,1)); 

figure(1), clf, h = gobjects(0); 
if any(named('-xla')), X_label = get_('-xla'); 
else X_label = {}; 
end
  
for ii = 1:size(X,2), px = X(:,ii);   
  if isa_cat_var(ii), px = px + randn(size(X(:,1)))/8; end
  
  subplot(size(X,2),1,ii)
  plot(px,Y,'o'), axis tight, tools.tidyPlot
  hold on, h(ii) = plot(px,0*Y,'+'); 
  if isa_cat_var(ii), set(gca,'XTick',min(X(:,ii)):max(X(:,ii))); end
  if ~isempty(X_label), 
    if iscell(X_label{ii}), set(gca,'XTickLabel',X_label{ii}); 
    else xlabel(X_label{ii}), 
    end
  end
end
subplot(size(X,2),1,1), title('equation')

if isempty(equations), 
  stats = dynamic_exploration(X,Y); 
  equations = stats.eqn;
  n_terms = stats.n_terms;
  labels = stats.models;
else 


  %% default: 6 terms in equation? 

  if any(named('-p0')), p0 = get_('-p0'); 
  else 
    p0 = [X X(:,1)*0+1] \ Y ;   
    if numel(p0) > 3, p0 = p0([4 3 1 1 2 2])' .* [1 1 1 0 1 0]; end
  end
  if max(n_terms) > numel(p0), p0(end+1 : max(n_terms)) = 0; end

  gof = zeros(size(equations));
  residuals = cell(size(equations)); 


  stats = struct; 
  stats.source = Y_label; 
  stats.models = labels;
  stats.par = cell(size(equations));
  stats.eqn = equations;
  stats.gof = gof; 

  if any(named('-aic')), stats.aic = gof; aic = gof; end

  fit_('optimset', f_opt); 

  for ii = 1:numel(equations) % fit each equation 

    % printInfo('[%s %2d] %s ',ty, ii, labels{ii})

    [p_est, gof(ii)] = fit_(equations{ii}, X, Y, p0); 

    Y_est = equations{ii}(X,p_est);
    stats.par{ii} = p_est; 

    assert(all(size(Y_est) == size(Y)))

    residuals{ii} = Y - Y_est;    

    set(h,'YData',Y_est); 
    title(sprintf('\\rm[%02d] %s : RMSE %0.2f', ii, labels{ii}, gof(ii)))
    ylabel(Y_label)

    if any(named('-aic'))
      [gof(ii,2),aic(ii,1),aic(ii,2)] = logLik(equations{ii},X,Y,p_est);
    end

    if do_PDF
         nom = sprintf('%sBpage-%02d-%03d.ps',ps_folder, ... 
                                              tt,ii);
         figs_to_ps(gcf,nom,'-move');
         pause(0.1), 
    else pause(0.1)
    end
  end % fit equtions loop

  %% Compute the grid of F tests

  stats.p = nan(numel(residuals)); 
  stats.F = nan(numel(residuals)); 
  stats.F_df = []; 
  stats.gof = gof; 
  if any(named('-aic')), stats.aic = aic; end

  [~,stats.residual_ks] = cellfun(@kstest,residuals);

  for ii = 2:numel(residuals)    
    for jj = 1:numel(residuals)
      [~,p,~,v] = vartest2(residuals{jj},residuals{ii},'tail','right');      
      stats.p(jj,ii) = p ;
      stats.F(jj,ii) = v.fstat;
      stats.F_df = [v.df1 v.df2];
    end
  end
  
end

nEq = numel(equations);

if any(named('-aic')), 
     gof = stats.aic(:,2); metric = 'AICc';
else gof = stats.gof(:,1); metric = 'rmse'; 
end

%%

if any(named('-fig')), figure(get_('-fig'))
else figure
end, clf, set(gcf,'Position', [180 150 1020 820])

img = log10(stats.p);
img(isnan(img)) = 0.2;
img((1:nEq) < (1:nEq)') = img((1:nEq) < (1:nEq)')/-min([img(:)*3; -8]);

imagesc(img), hold on

if nEq > 32
  set(gca,'XAxisLocation','top','XTickLabelRotation',-45)
else
  set(gca,'XTick',1:nEq,'XTickLabel',labels, ...
          'YTick',1:nEq,'YTickLabel',labels, ...
          'XAxisLocation','top','XTickLabelRotation',-45)
end
axis image, colorbar, colormap(bone)
caxis([-6 0])

C = lines(7); W = @(i,a) (C(i,:)+(1/a)) / (1 + 1/a);
color = interp1([min(gof); median(gof); max(gof)], C([5 3 2],:), gof);
scatter(1:nEq, 1:nEq, [], color, 'o','filled')

u = pi*linspace(0,1,31)/2; 

xy_edge = [0 0; [1 0] + 0.3*[(sin(u')-1) (1-cos(u'))]; 1 1];
xy_edge = @(i,j,d) xy_edge(:,d)*(i-j) + j; 

% plot(xy_edge(3,5,1), xy_edge(3,5,2), 'g-','LineWidth',3)

%%  

fprintf('\n===========================\n%s:\n', Y_label)

alpha = 0.05 ; 

this_picked = false(size(n_terms)); 
this_picked(n_terms == 1) = 1; 

tps_idx = 1; 

for nn = 1:(max(n_terms))

  up_1 = (n_terms == nn+1); 
  % assert(any(up_1));

  if ~any(this_picked(nn)), 
       n_eff = max(n_terms(this_picked)); 
  else n_eff = nn; 
  end

  for ii = 1:numel(equations)

    if n_terms(ii) ~= n_eff, continue, end
    if ~this_picked(ii), continue, end

    if ~any(named('-q'))    
      fprintf('\n[%02d] %9.2f %s\n', ii, gof(ii), labels{ii})
    end
    iidx = find(stats.p(ii,:) < alpha & up_1');         
    % if nn == 1, iidx = find(up_1) ; end

    text(ii+0.5,ii,sprintf('%s = %0.2f',metric,gof(ii)),'Color', color(ii,:))

    if ~any(up_1), continue, end

    if ii == tps_idx && any(iidx) % update 'best model'
      [~,seq] = nanmin(stats.p(ii,iidx));
      tps_idx = iidx(seq); 
    end

    [~,seq] = sort(gof(iidx),'ascend');

    if ~isempty(seq)
      for ss = seq', jj = iidx(ss);
        p = stats.p(ii,jj); 
        if ~any(named('-q'))    
          fprintf(' --> p=%6.5f [%02d] %s\n', p, jj, labels{jj})
        end
        plot(xy_edge(ii,jj,1), xy_edge(ii,jj,2), ...
                  'Color', [C(2,:) min(1,0.01/p)], ...
                  'LineWidth',5*min(stats.p(ii,iidx))/p)  
        text(jj-0.15,jj+0.5,sprintf('p = %0.4f',p),'Color', W(2,min(1,0.01/p)),'Fontsize',8,'Horiz','right')
      end
    end

    if ~any(named('-q'))    
      fprintf(' (%d models n/s)\n', sum(up_1) - length(iidx))
    end

    this_picked(iidx) = true; 
  end
end

scatter(1:nEq, 1:nEq, 100, color, 'o','filled')   
title(Y_label,'FontSize',20)
% set(gca,'Position',get(gca,'Position')-[0 0.03 0 0])


if do_PDF
     nom = sprintf('%sApage-%02d.ps',ps_folder,tt);
     figs_to_ps(gcf,nom,'-move');
     pause(0.1), 
else pause(0.1)
end

%%

stats.possible_model = this_picked;  
stats.selected_model = tps_idx;

if any(named('-aic'))  
  stats.F_test_selected_model = tps_idx;
  [~,stats.selected_model] = min(stats.aic(:,2));
end

ii = tps_idx; 
fprintf('\nBest model:\n[%02d] %s\n', ii, labels{ii})

fprintf('p > %0.4f\n',min(stats.p(ii,n_terms > n_terms(ii))))

if ~any(named('-no-boot')), stats = bootstrap_parameters(X, Y, stats); end

return


function [equations, labels, n_terms] = setup_statistical_models

named = evalin('caller','named');
get_  = evalin('caller','get_');


re_label = {'@(x,p)',''}; 
isa_term = {}; 

if any(named('-eq')), equations = get_('-eq');
elseif any(named('-peri')),  

  do_v2 = any(named('-standard-form')) || any(named('-v2'));   
  equations = stats_models_peri_thickness(do_v2); 

  re_label = {'@(x,p)',''; 'p(1)','c'; 'p_xy(x,p)','c_{xy}'; '.*','*'; ...
                           'p(3)','a_1'; 'p_xy(x,p(3:4))','a_{1,xy}'; ...
                           'p(5)','a_2'; 'p_xy(x,p(5:6))','a_{2,xy}'; ...
                           'p(7)','a_{12}'; 'p_xy(x,p(7:8))','a_{12,xy}'; ...
                           'x(:,1)','\bfPW\rm'; 'x(:,2)','\bfFD\rm' };

  isa_term = {'y','a','c','7'}; 
else

  % this approach falls over for large nX
  
  equations = {}; labels = {}; n_terms = []; 
  return
%     
%     
%   nX = evalin('caller','size(X,2)');
%   
%   if nX > 5, 
%   end
%   
%   if any(named('-all-int')), error all_intersections, end
%    
%   equations = { @(x,p) p(1) }; 
% 
%   nTerm = 1; 
%   while numel(equations) < 1e3 
%     nex_layer = nchoosek(1:nX,nTerm);
%     
%     for ii = 1:size(nex_layer,1)
%       
%       coe = [(1:size(nex_layer,2)); nex_layer(ii,:)];
%       eqn = sprintf(' + p(%d)*x(:,%d)',coe(:));      
%       equations{end+1} = str2func(sprintf('@(x,p) p(1)%s',eqn));
%       
%       if nTerm  == 1, continue, end
%       ix_terms =  nchoosek(nex_layer(ii,:),2);
%       coe = [(1:size(nex_layer,2)); nex_layer(ii,:)];
% 
%       
%       
%       eqn2 = 
%       
%     end
%     
%     if nTerm > 1, 
%       ix_terms =  nchoosek(1:nX,2);
% 
%     
%     
%     for ii = 1:size(nex_layer,1)
%       
%       coe = [(1:size(nex_layer,2)); nex_layer(ii,:)];
%      
%       
%       eqn = sprintf(' + p(%d)*x(:,%d)',coe(:));      
%       equations{end+1} = str2func(sprintf('@(x,p) p(1)%s',eqn));
%     end
% 
%     if nTerm == nX, break, end
%     nTerm = nTerm + 1; 
%   end
%   
%   
%   
%   
  
% 
%   error set_up_equations_from_X_Y
end

%% From the block of candidate equations get display labels and n_terms: 
if any(named('-lbl')), labels = get_('-lbl'); 
else  
  if any(named('-re')), re_label = get_('-re'); end
  
  labels = cellfun(@func2str,equations,'unif',0); 
  for ii = 1:size(re_label,1), labels = strrep(labels, re_label{ii,:}); end
  % labels{1} = 'c'; 
end

if any(named('-te')), isa_term = get_('-te'); end

n_terms = @(eq) sum(cellfun(@(t) sum(eq == t), isa_term)); 
n_terms = cellfun(n_terms, labels);

[n_terms,seq] = sort(n_terms,'ascend');
labels = labels(seq);
equations = equations(seq); 


% for large # of columns, equations explode if determined in advance
function stats = dynamic_exploration(X,Y) % in-context function

named = evalin('caller','named');
get_  = evalin('caller','get_');

%   equations = { @(x,p) p(1) }; 
% 
%   nTerm = 1; 
%   while numel(equations) < 1e3 
%     nex_layer = nchoosek(1:nX,nTerm);
%     
%     for ii = 1:size(nex_layer,1)
%       
%       coe = [(1:size(nex_layer,2)); nex_layer(ii,:)];
%       eqn = sprintf(' + p(%d)*x(:,%d)',coe(:));      
%       equations{end+1} = str2func(sprintf('@(x,p) p(1)%s',eqn));
%     end
% 
%     nTerm = nTerm + 1; 
%   end

equations = { @(x,p) p(1) };

nX = size(X,2); 

for ii = 1:nX,   
  term = sprintf(' + p(%d)*x(:,%d)',2,ii);      
  equations{end+1,1} = str2func(sprintf('@(x,p) p(1)%s',term));
end

do_PDF = any(named('-pdf')); 

active_pre = (0:nX)' == 0;
active_post = (0:nX)' > 0;

max_order = 2*nX; 
if any(named('-max-order')), max_order = get_('-max-order'); end

gof = zeros(size(equations));
residuals = cell(size(equations)); 

stats = struct; 
stats.source = evalin('caller','Y_label'); 
stats.models = cellfun(@func2str,equations,'unif',0);
stats.par = {};
stats.eqn = equations;
stats.gof = gof; 

stats.p = []; 
stats.F = []; 

if any(named('-aic')), stats.aic = gof; aic = gof; end

fit_('optimset', evalin('caller','f_opt')); 

h = evalin('caller','h'); 


if any(named('-cv'))  
  
  if any(named('-kf')),   
       chunks = cvpartition(size(X,1),'KFold',get_('-kf'));
  else chunks = cvpartition(size(X,1),'KFold',10);  
  end
end

while true 

  for ii = 1:numel(equations) % fit each equation 

    if numel(residuals) >= ii && ~isempty(residuals{ii}), continue, end
    
    eq = func2str(equations{ii}); 
    x_in = arrayfun(@(x) contains(eq,sprintf('x(:,%d)',x)), 1:nX); 
    n_p = max(str2double(regexp(eq,'(?<=p\()\d+','match')));
    
    p0 = [X(:,1)*0+1 X(:,x_in)] \ Y; 
    p0(end+1 : n_p) = 0;

    % printInfo('[%s %2d] %s ',ty, ii, labels{ii})
    
    if any(named('-cv'))
      
      Y_est = 0*Y; 
      p_cv = []; 
      
      for cc = 1:chunks.NumTestSets
      
        [p_cv(cc,:), ~] = fit_(equations{ii},X(chunks.training(cc),:),Y(chunks.training(cc)),p0,false); %#ok<AGROW>
        Y_est(chunks.test(cc)) = equations{ii}(X(chunks.test(cc),:),p_cv(cc,:));      
      end     
      gof(ii,1) = sqrt(nanmean((Y - Y_est).^2));
      stats.par{ii,1} = nanmean(p_cv,1);
      
    else
      [p_est, gof(ii,1)] = fit_(equations{ii}, X, Y, p0, false); 
      stats.par{ii,1} = p_est; 

      Y_est = equations{ii}(X,p_est);     
      if numel(Y_est) == 1, Y_est = Y_est * ones(size(Y)); end
    end
      
    assert(all(size(Y_est) == size(Y)) || numel(Y_est) == 1)
    residuals{ii} = Y - Y_est;    

    if ishandle(h), set(h,'YData',Y_est); end    
    title(sprintf('\\rm[%02d] %s : RMSE %0.2f', ii, eq, gof(ii,1)))
    ylabel(stats.source)

    if any(named('-aic'))
      [gof(ii,2),aic(ii,1),aic(ii,2)] = logLik(equations{ii},X,Y,p_est);
    end

    if do_PDF
         nom = sprintf('%sBpage-%02d-%03d.ps',ps_folder, ... 
                                              tt,ii);
         figs_to_ps(gcf,nom,'-move');
         pause(0.1), 
    else pause(0.1)
    end
  end % fit equtions loop

  for ii = 2:numel(residuals)    
    for jj = 1:numel(residuals)
      if all(size(stats.F) >= [jj ii]) && stats.F(jj,ii) > 0, continue, end
      [~,p,~,v] = vartest2(residuals{jj},residuals{ii},'tail','right');      
      stats.p(jj,ii) = p ;
      stats.F(jj,ii) = v.fstat;
    end
  end
  
  %% Based on F test add new equations to set
  
  if n_p == max_order, break, end
  
  seeds = active_post; 
  seeds(seeds) = any(stats.p(active_pre,active_post) < 0.05,1);
  
  if sum(seeds) > nX/2 % limit growth
    seeds(seeds) = gof(seeds,1) > min(maxk(gof(seeds,1),ceil(nX/2)));
  end  
  
  for ii = find(seeds)'
    eq = func2str(equations{ii}); 
    x_in = arrayfun(@(x) contains(eq,sprintf('x(:,%d)',x)), 1:nX); 
    n_p = max(str2double(regexp(eq,'(?<=p\()\d+','match')));
    
    for jj = 1:nX,
      if x_in(jj), continue, end
      term = sprintf('%s + p(%d)*x(:,%d)',eq,n_p+1,jj);
      
      term = strtrim(strsplit(term,'+')); % put in x() order
      [~,seq] = sort(cellfun(@fliplr,term(2:end),'unif',0));
      term = sprintf('%s + ',term{[1 seq+1]}); 
      term(end-2:end) = ''; % trim final + ..
      
      if any(strcmp(stats.models,term)), continue, end % already added      
      equations{end+1,1} = str2func(term); %#ok<AGROW>
      stats.models{end+1,1} = term;
    end
    
    if sum(x_in) == 1, interactions = [];
    else interactions = nchoosek(find(x_in),2);      
    end
    for jj = 1:size(interactions,1)
      
      term = sprintf('%s + p(%d)*x(:,%d).*x(:,%d)',eq,n_p+1,interactions(jj,:));
      
      term = strtrim(strsplit(term,'+')); % put in x() order
      [~,seq] = sort(cellfun(@fliplr,term(2:end),'unif',0));
      term = sprintf('%s + ',term{[1 seq+1]}); 
      term(end-2:end) = ''; % trim final + ..
      
      if any(strcmp(stats.models,term)), continue, end % already added      
      equations{end+1,1} = str2func(term); %#ok<AGROW>
      stats.models{end+1,1} = term;
    end
  end
  
  active_pre = seeds; 
  active_post(:) = false; active_post(end+1 : numel(equations)) = true; 
  if numel(active_pre) == numel(active_post), break, end
  
end

[~,stats.residual_ks] = cellfun(@kstest,residuals);

stats.eqn = equations;


stats.gof = gof; 
if any(named('-aic')), stats.aic = aic; end

stats.n_terms = cellfun(@(eq) max(str2double(regexp(eq,'(?<=p\()\d+','match'))), stats.models);

labels = stats.models; % make nicer formatting
labels = regexprep(labels,'\((:,)?(\d+)\)','_{$2}'); 
labels = regexprep(labels,'\.?*',' ');
labels = strrep(labels,'@(x,p)','');
stats.models = labels; 

return


function equations = stats_models_peri_thickness(do_v2)



p_xy = @(x,p) p(1) + p(2).*x(:,3); 

equations = { @(x,p) p(1) + 0*x(:,1)  % constant (+/- effect of fasc location)
              @(x,p) p_xy(x,p)
              @(x,p) p(1) + p(3) .* x(:,1)               % factor 1 : peri thickness
              @(x,p) p_xy(x,p) + p(3) .* x(:,1)
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)) .* x(:,1)              
              @(x,p) p(1) + p(5) .* x(:,2)               % factor 2 : fasc diam
              @(x,p) p_xy(x,p) + p(5) .* x(:,2)
              @(x,p) p_xy(x,p) + p_xy(x,p(5:6)) .* x(:,2)              
              @(x,p) p(1) + p(3) .* x(:,1) + p(5) .* x(:,2) % both factors 1 and 2 (additive)
              @(x,p) p_xy(x,p) + p(3) .* x(:,1) + p(5) .* x(:,2)
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)) .* x(:,1) + p_xy(x,p(5:6)) .* x(:,2)
              @(x,p) p(1) + (p(3) + p(5).*x(:,2) ) .* x(:,1);  % both factors 1 and 2 (multiplicative, 1-dominant)
              @(x,p) p_xy(x,p) + (p(3) + p(5).*x(:,2) ) .* x(:,1)
              @(x,p) p_xy(x,p) + (p_xy(x,p(3:4)) + p(5).*x(:,2) ) .* x(:,1)
              @(x,p) p_xy(x,p) + (p(3) + p_xy(x,p(5:6)).*x(:,2) ) .* x(:,1)
              @(x,p) p_xy(x,p) + (p_xy(x,p(3:4)) + ...
                                  p_xy(x,p(5:6)).*x(:,2) ) .* x(:,1)
              @(x,p) p(1) + (p(5) + p(3).*x(:,1) ) .* x(:,2)  % both factors 1 and 2 (multiplicative, 2-dominant)
              @(x,p) p_xy(x,p) + (p(5) + p(3).*x(:,1) ) .* x(:,2)
              @(x,p) p_xy(x,p) + (p_xy(x,p(5:6)) + p(3).*x(:,1) ) .* x(:,2)
              @(x,p) p_xy(x,p) + (p(5) + p_xy(x,p(3:4)).*x(:,1) ) .* x(:,2)
              @(x,p) p_xy(x,p) + (p_xy(x,p(5:6)) + ...
                                  p_xy(x,p(3:4)).*x(:,1) ) .* x(:,2) };

                                
if do_v2
                          
equations = { @(x,p) p(1) + 0*x(:,1)  % constant (+/- effect of fasc location)
              @(x,p) p_xy(x,p)          % factor 3 : XY seperation
              @(x,p) p(1) + p(3)*x(:,1) % factor 1 : peri thickness
              @(x,p) p(1) + p(5)*x(:,2) % factor 2 : fasc diam
              @(x,p) p_xy(x,p) + p(3)*x(:,1) % factor 1,3 : peri thickness
              @(x,p) p_xy(x,p) + p(5)*x(:,2) % factor 2,3 : fasc diam
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)).*x(:,1) % factor 3,1,13
              @(x,p) p_xy(x,p) + p_xy(x,p(5:6)).*x(:,2) % factor 3,2,23
              @(x,p) p(1) + p(3) .* x(:,1) + p(5) .* x(:,2) % both factors 1 and 2 (additive), no factor 3
              @(x,p) p_xy(x,p) + p(3) .* x(:,1) + p(5) .* x(:,2) % all factors 1 - 3 (additive), no interactions
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)).*x(:,1)+p_xy(x,p(5:6)).*x(:,2) % factor 1-3,13,23
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)).*x(:,1)+p(5).*x(:,2)+p(7).*x(:,1).*x(:,2) % all factor 1 interactions
              @(x,p) p_xy(x,p) + p(3).*x(:,1)+p_xy(x,p(5:6)).*x(:,2)+p(7).*x(:,1).*x(:,2) % all factor 2 interactions
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)).*x(:,1)+p_xy(x,p(5:6)).*x(:,2) % all factor 3 interactions
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)).*x(:,1)+p_xy(x,p(5:6)).*x(:,2) + p(7).*x(:,1).*x(:,2) % all 2-term interactions               
              @(x,p) p_xy(x,p) + p_xy(x,p(3:4)).*x(:,1)+p_xy(x,p(5:6)).*x(:,2) + p_xy(x,p(7:8)).*x(:,1).*x(:,2) }; % all interactions 
              
end

function [pFit, gof] = fit_(eq, X, Y, pInit, do_several)

persistent f_opt once
if ischar(eq), f_opt = X; once = true; return, end
if nargin < 5, do_several = true; end

lsq_ = @(p,f) nanmean((Y - f(X,real(p))).^2);  
root_ = @(p,idx) p.*(~idx) + sqrt(p).*idx; 

% printInfo('[%s %2d] %s ',ty, ii, labels{ii})

% Fit the equation from a few different starting points 
[pFit(1,:),sel(1)] = fminsearch(@(p)lsq_(p,eq), pInit, f_opt);

if do_several
  if numel(pInit) == 6
    [pFit(2,:),sel(2)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 1 1 0 0]), f_opt); 
    [pFit(3,:),sel(3)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 1 1 1 1]), f_opt);     
    [pFit(4,:),sel(4)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 0 0 1 1]), f_opt); 
  elseif numel(pInit) == 8
    [pFit(2,:),sel(2)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 1 1 0 0 0 0]), f_opt); 
    [pFit(3,:),sel(3)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 1 1 1 1 0 0]), f_opt);     
    [pFit(4,:),sel(4)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 0 0 1 1 0 0]), f_opt); 
    [pFit(5,:),sel(5)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 1 1 0 0 1 1]), f_opt); 
    [pFit(6,:),sel(6)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 1 1 1 1 1 1]), f_opt);     
    [pFit(7,:),sel(7)] = fminsearch(@(p)lsq_(p,eq), root_(pInit,[0 0 0 0 1 1 1 1]), f_opt);
  elseif once, warning('FtestSeq:augment:unclear', ...
                       'Unknown # of parameters %d', numel(pInit));
         once = false;
  end
end

[~,sel] = min(sel); 
pFit = real(pFit(sel,:));
  
gof = sqrt(lsq_(pFit, eq)); 

function [lnl,aic,aicc] = logLik(eq,X,Y,par)

if isnumeric(eq), rem = X; par = Y;
else rem = eq(X,par) - Y; 
end

s = std(rem);
n = numel(rem);

lnl = -numel(rem) * log(s * sqrt(2*pi)) - sum( (rem / s).^2 / 2);

if nargout < 2, return, end

if numel(par) > 1
  % which variables are even relevent to this model?
  if any(isnan(eq(X,par))), check = @all; else check = @any; end
  for ii = 1:numel(par)
    p_chk = par;
    p_chk(ii) = nan;
    relevent(ii) = check(isnan(eq(X,p_chk)));
  end  
  n_par = sum(relevent) + 1; % because sigma counts as a param
else n_par = par + 1; 
end

aic = -2*lnl + 2*n_par; 
aicc = aic + 2*n_par*(n_par+1)/(n-n_par-1); 

function S = bootstrap_parameters(X,Y,S) 

named = evalin('caller','named');
get_  = evalin('caller','get_');


if any(named('-nb')), nB = get_('-nb'); 
else nB = 1e3; 
end

eq = S.eqn{S.selected_model};
p0 = S.par{S.selected_model};

pX = zeros(nB,1) * reshape(p0,1,[]); 

fprintf('Running bootstrap (%d iterates) for Y=%s \n', nB, S.models{S.selected_model})

for rep = 1:nB  
  idx = ceil(rand(size(Y))*numel(Y));  
  pX(rep,:) = fit_(eq,X(idx,:),Y(idx,:),p0.*(randn(size(p0))/2+1), false);  
end

S.bootstrap.model = S.selected_model;
S.bootstrap.n_replicates = nB; 
S.bootstrap.variables = true(size(p0)); 
S.bootstrap.mean = mean(pX);
S.bootstrap.std = std(pX);
[~,~,S.bootstrap.ci] = ttest(pX);

if any(named('-boot-qq')), qq = get_('-boot-qq'); else qq = [1:3]/4; end
S.bootstrap.QQ = reshape(qq,1,[]);
S.bootstrap.IQR = quantile(pX,S.bootstrap.QQ);

% which variables are even relevent to this model?
if any(isnan(eq(X,p0))), check = @all; else check = @any; end

for ii = 1:numel(p0)
  
  p_chk = p0;
  p_chk(ii) = nan;
  
  S.bootstrap.variables(ii) = check(isnan(eq(X,p_chk)));
end






















