
function make_axon_population(axons_file, ...
                              nerve_xml, ...
                              nerve_script, ... 
                              n_axon_unit, ...
                              n_myelinated, ...
                              n_unmyelinated, ...
                              aff_eff_ratio, varargin ) 

if nargin < 1 || isempty(axons_file)
  if nargin == 0, 
     n_myelinated = 70; n_unmyelinated = 300; n_axon_unit = 'count';
     fprintf('*** running in demo mode\n')
  end
  
  fprintf('Arg 1 not set, using default axons.mat file (rat-cervical-vagus)\n')
  axons_file = 'rat-cervical-vagus'; 
end
if nargin < 2
  fprintf('Arg 2 not set, using default nerve.xml file\n')
  nerve_xml = './input/demo/nerve.xml'; 
end
if nargin < 3, nerve_script = './input/demo/nerve-script.json';
end

named = @(v) strncmpi(v,varargin,length(v)); 
has_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if any(named('-q')), printf = @(s,varargin) 0; 
else printf = @(s,varargin) fprintf([s '\n'],varargin{:}); 
end

if true
    disp('====================================================')
    fprintf('Running models.axon_population %s\n', datestr(now))
    fprintf('{1} = %s\n{2} = %s\n{3} = %s\n',axons_file, nerve_xml, nerve_script)
    disp('====================================================')
end

t = tic; 

tools.file('root',pwd); % set 'root' to this folder
if ~exist(axons_file,'file'), axons_file = ['./input/' axons_file '.mat']; end
printf('Loading %s', axons_file);
ax = load(axons_file); 

if isdeployed, disp('Progress: 5%'), end

%% Figure out nerve anatomy 
if ~isempty(nerve_script)
 
  printf('Loading %s', nerve_script);
  if has_ext_(nerve_script,'.mat'), n = load(nerve_script); 
  elseif has_ext_(nerve_script,'.json'), n = tools.parse_json(nerve_script);
  else error('%s: unknown extension, expected .mat or .json', nerve_script)
  end      

  if isfield(n,'nerve'), opts.nerve = n.nerve; else opts.nerve = n; end
  if isfield(n,'mesh'), opts.mesh = n.mesh; end
end

if ~isempty(nerve_xml), opts.nerve.file = nerve_xml;
elseif ~isfield(opts,'nerve'), opts.nerve.file = nerve_xml;
end


mesh.insert_gmsh_fascicles('-setup',opts);
nerve = mesh.insert_gmsh_fascicles('-info','-check-units');

if isfield(nerve,'splines') && numel(nerve.splines) > 1
    nsn = sprintf(', %s',nerve.splines.type); 
    fprintf('selecting %s from {%s} in %s\n', nerve.splines(1).type, ...
                                              nsn(3:end), opts.nerve.file)
    nerve.splines = nerve.splines(1); 
end

if isdeployed, disp('Progress: 15%'), end

output_file = tools.file('get','./output/axon-population (%d).mat','next');
arguments = {'-anat',nerve,'-out', output_file}; 

while 1
  if ~exist('n_axon_unit','var'), break, end
  if nargin < 7, aff_eff_ratio = 4; end % Foley and DuBois, 1937
  if ischar(aff_eff_ratio),  aff_eff_ratio = str2double(aff_eff_ratio);   end
  if ischar(n_myelinated),   n_myelinated  = str2double(n_myelinated);    end
  if ischar(n_unmyelinated), n_unmyelinated = str2double(n_unmyelinated); end
  
  eff_frac = 1/(1 + aff_eff_ratio); 
  aff_frac = 1 - eff_frac;
  n_axons  = [ [aff_frac eff_frac] * n_myelinated ...
               [aff_frac eff_frac] * n_unmyelinated ];
  switch lower(n_axon_unit)
    case 'per_mm2'
      try
        nF = size(nerve.splines.outline,3);
      catch C, disp(nerve), disp(nerve.splines), 
        error('Crashed at nerve.splines.outline')
      end
      f_area_mm2 = arrayfun(@(f) polyarea(nerve.splines.outline(:,1,f), ...
                                          nerve.splines.outline(:,2,f)), ...
                                   1:nF);
      n_axons = round(n_axons * f_area_mm2);
    case 'count',  n_axons = round(n_axons);
    case 'ignore', n_axons = cellfun(@numel,{ax.axon_populations.fibre_diam});
    otherwise warning('unknown n_axon_unit value %s', n_axon_unit)
        n_axons = []; 
  end
  if ~isempty(n_axons)
    arguments = [arguments {'-counts', n_axons }]; %#ok<AGROW>
    printf('Axon Counts = [%s ]', sprintf(' %d',arguments{end}));
  end
  break
end


axon_types = unique({ax.axon_populations.axon_model});
currents = struct; 
for ty = axon_types
  if isfield(ax,ty{1}), currents.(ty{1}) = ax.(ty{1}); end
end
arguments = [arguments {'-simulations',currents}];

if any(named('-:')), arguments = [arguments ...
                                  varargin(find(named('-:'))+1:end)]; 
end
axon_population('-pregenerated',ax.axon_populations,arguments{:});

if isdeployed, disp('Progress: 90%'), end

pause(0.05)
if ishandle(2), print(2,strrep(output_file,'.mat','.svg'),'-dsvg'), end
pause(0.05), close all

if isdeployed, disp('Progress: 100%'), end

toc(t)


return

function axon_population(varargin)
% make_axon_population takes a set of fascicle contours and a specification
% for the axon populations, and populates the fascicles with axons with
% varying diameters and g-ratios. 
% 
% We have five data-sets: four scraped from published papers and
%  one which was provided by SPARC. Data-sets scraped from published papers
%  rely on reconstructing the underlying joint distribution from the 
%  published marginal histograms.
% 
% NOTE the remainder of this help code is moderately out-of-date
% 
% make_axon_population -sol uses N Soltanpour and RM Santer, 
%   "Preservation of the cervical vagus nerve in aged rats," (1996)
% make_axon_population -pow uses JC Prechtl and TL Powley, 
%   "The fiber composition of the abdominal vagus of the rat," (1990)
% make_axon_population -pn uses N Biscola, L Haveton (unpublished)
% 
% the default contours and axon class counts are drawn from figure 2 and 
%   table 6 of CE Hulsebosch and R Coggeshall, "An analysis of the axon 
%   populations in the nerves to the pelvic viscera in the rat" (1982)
% 
% alternate set of contours can be supplied by calling 
%   nerve = mesh.read_dat_file('filename.xml');
%   models.make_axon_population( ..., '-in', nerve );
% 
% alternate target nerve counts can be supplied by calling 
%   make_axon_population( ..., '-class', table ); 
%   where table is a 2x2 table: [A-type, C-type afferents (sensory);
%                                A-type, C-type efferents (motor)  ]
%   The default (Hulsebosch) is: [  475, 1200; 
%                                   395, 2800  ]
% 
% By default, two figures are generated. These can be suppressed by calling
%   make_axon_population ... -no-fig
% 
% Version 0.2, refactored to enable code-free addition of new datasets 
% Version 0.2, refactored to merge approaches
% Calvin Eiber 16-April-2020 ceiber@unimelb.edu.au

varargin = tools.opts_to_args(varargin,'axons');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

D = struct; 
assert(any(named('-pregen')) && any(named('-simulation')), ...
       'This iteration of this code works with pre-generated axon distributions')

axons = get_('-pregen');

if any(named('-class')), tab = get_('-class');
elseif any(named('-counts')), tab = get_('-counts'); 
else tab = cellfun(@numel,{axons.fibre_diam});
end

t = struct;
t.Type = repmat({'unknown type'},[numel(tab) 1]);
t.Count = tab(:);
t.Model = repmat({''},[numel(tab) 1]);
if numel(tab) >= 1, t.Type{1} = 'Myelinated Afferent'; 
                    t.Model{1} = 'Gaines';              end
if numel(tab) >= 2, t.Type{2} = 'Myelinated Efferent'; 
                    t.Model{2} = 'MRG';                 end
if numel(tab) >= 3, t.Type{3} = 'Unmyelinated Afferent'; 
                    t.Model{3} = 'Sundt';               end
if numel(tab) >= 4, t.Type{4} = 'Unmyelinated Efferent'; 
                    t.Model{4} = 'Sundt';               end
tab = struct2table(t);

nerve = get_anatomy(get_,named);

for ty = 1:numel(axons) % up or down sample 
  
  n = tab.Count(ty); 
  nP = numel(axons(ty).fibre_diam);

  if numel(axons(ty).fibre_diam) >= n 
    [~,subset] = sort(rand(nP,1)); 
    for f = fieldnames(axons)'
      if size(axons(ty).(f{1}),1) ~= nP, continue, end
      axons(ty).(f{1}) = axons(ty).(f{1})(subset(1:n),:);
    end
  else
    superset = ceil(nP*rand(1,n-nP));
    for f = fieldnames(axons)'
      if size(axons(ty).(f{1}),1) ~= nP, continue, end
      axons(ty).(f{1}) = axons(ty).(f{1})([1:end superset],:);
    end
  end
end


temp = struct; 

temp.fibre_diam = [axons(1).fibre_diam; axons(2).fibre_diam];
temp.unmyelinated_diam = [axons(3).fibre_diam; axons(4).fibre_diam];

temp = arrange_axons(D, nerve, temp, tab); 

axons(1).fascicle = temp.fibre_fascicle(temp.fibre_afferent);
axons(2).fascicle = temp.fibre_fascicle(~temp.fibre_afferent);
axons(3).fascicle = temp.unmyelinated_fascicle(temp.unmyelinated_fascicle);
axons(4).fascicle = temp.unmyelinated_fascicle(~temp.unmyelinated_fascicle);

axons(1).axon_xy = temp.fibre_xy(temp.fibre_afferent,:);
axons(2).axon_xy = temp.fibre_xy(~temp.fibre_afferent,:);
axons(3).axon_xy = temp.unmyelinated_xy(temp.unmyelinated_fascicle,:);
axons(4).axon_xy = temp.unmyelinated_xy(~temp.unmyelinated_fascicle,:);

%%
if ~any(named('-no-fig'))
    %%
    nF = max(cat(1,axons.fascicle));
    color = cat(1,axons.color);

    figure(2), clf, hold on, % C = lines(nF);
    set(gcf,'Position',get(gcf,'Position') .* [1 .7 1.3 1])
    if iscell(nerve.outline)
        fill(nerve.outline{end}(:,1),nerve.outline{end}(:,2),[.9 .9 .9],'EdgeColor','none')
    end  
    for ff = 1:nF
      fill(nerve.fascicles(:,1,ff),nerve.fascicles(:,2,ff),'w','EdgeColor',[.3 .3 .3],'LineWidth',1.2)    

      sel = (f_id == ff) & ~is_myelin & is_afferent;
      plot(xy(sel,1),xy(sel,2),'.','Color',color(1,:),'MarkerSize',10) 

      sel = (f_id == ff) & ~is_myelin & ~is_afferent;
      plot(xy(sel,1),xy(sel,2),'.','Color',color(2,:),'MarkerSize',10) 

      sel = (f_id == ff) & is_myelin & is_afferent;
      plot(xy(sel,1),xy(sel,2),'o','MarkerFaceColor',(color(3,:)+1.5)/2.5, ...
                  'MarkerSize',5,'Color',color(3,:),'LineWidth',1.2)

      sel = (f_id == ff) & is_myelin & ~is_afferent;
      plot(xy(sel,1),xy(sel,2),'o','MarkerFaceColor',(color(4,:)+1.5)/2.5, ...
                  'MarkerSize',5,'Color',color(4,:),'LineWidth',1.2)
    end
    axis image, tools.tidyPlot

    l_str = {'Fascicle outline','unmyelinated afferent', ...
             'unmylinated efferent','myelinated afferent', ...
                                    'myelinated efferent'};
    for tt = height(tab):-1:1    
    if tab.Count(tt) == 0, l_str(tt+1) = []; end
    end

    if iscell(nerve.outline), l_str = [{'Endoneurium'} l_str]; end
    legend(l_str{:},'location','best')
end
%%
  
if any(named('-out')), output = get_('-out'); 
elseif any(named('-file')), output = get_('-file'); 
else output = tools.file('get','out~/axons (%d).mat','next');
end
if ~any(ismember(output,'/\'))
    output = ['sub~/axons/axons (' output ').mat'];
end

if any(output == '~'), output = tools.file(output); end

axon_models = get_('-sim'); % data pass-through

if any(named('-csv')), export_tables(axons,sam,nerve,'csv',output)  
elseif any(named('-xls')), export_tables(axons,sam,nerve,'xlsx',output)
else  
  if ~isfolder(fileparts(output)), mkdir(fileparts(output)); end
  while exist(output,'file'), output = strrep(output,'.mat','_NEW.mat'); end
  fprintf('Saving %s\n',tools.file('T',output))
  save(output,'axons','nerve','axon_models');
end

     %%
return

function pop = arrange_axons(D, nerve, pop, axon_types)

%% Generate random axon arrangements 
nA = sum(axon_types.Count); 

nerve.f_bounds      = min(permute(min(nerve.fascicles,[],1), [2 3 1]),[],2)';
nerve.f_bounds(2,:) = max(permute(max(nerve.fascicles,[],1), [2 3 1]),[],2)';

% Pass 1 - random scatter to get baseline fill ratio
xy = rand(nA,2) .* [nerve.f_bounds(2) - nerve.f_bounds(1) ...
                    nerve.f_bounds(4) - nerve.f_bounds(3)] + ...
                   [nerve.f_bounds(1)   nerve.f_bounds(3)];

nF = size(nerve.fascicles,3);
f_id = arrayfun(@(n) in_loop(xy,nerve.fascicles(:,:,n)), 1:nF,'Unif',0);
f_id = [f_id{:}]; 

in_a_bv = false(size(f_id(:,1))); % In a blood vessel?
if iscell(nerve.outline)
  for bb = 1:size(nerve.outline{2},3)
    in_a_bv = in_a_bv | in_loop(xy,nerve.outline{2}(:,:,bb));
  end
end

% Pass 2 - overgenerate to get 1.1x as many as we should need 
n_ratio = mean(any(f_id,2) & ~in_a_bv); 

while true % might need multiple passes 
  xy = rand(ceil(1.1*nA/n_ratio),2) .* [nerve.f_bounds(2) - nerve.f_bounds(1) ...
                                         nerve.f_bounds(4) - nerve.f_bounds(3)] + ...
                                        [nerve.f_bounds(1)   nerve.f_bounds(3)];
  for lloyd_relax = 1:3 % make more uniform 
     [cp,vp] = voronoin(xy); 
      cp(1,:) = NaN;    
      xy = [cellfun(@(p) nanmedian(cp(p,1)), vp) ...
            cellfun(@(p) nanmedian(cp(p,2)), vp)];

      if isdeployed, fprintf('Progress: %0.1f%%\n', 15 + (lloyd_relax/3)*70), end
      % h.XData = xy(:,1); h.YData = xy(:,2);
  end

  f_id = arrayfun(@(n) in_loop(xy,nerve.fascicles(:,:,n)), 1:nF,'Unif',0);
  f_id = [f_id{:}];  ok = any(f_id,2);

  in_a_bv = false(size(f_id(:,1))); % repeat "in a blood vessel"
  if iscell(nerve.outline)
    for bb = 1:size(nerve.outline{2},3)
      in_a_bv = in_a_bv | in_loop(xy,nerve.outline{2}(:,:,bb));
    end
  end
  ok(in_a_bv) = 0; 

  vp(~ok,:) = []; xy(~ok,:) = []; f_id(~ok,:) = []; % Remove extrafascicular
  [f_id,~] = find(f_id'); % convert to #
  f_id = reshape(f_id,[],1); % column vector

  if size(xy,1) < nA, % decrease n_ratio and try again
    if n_ratio < 0.01,
      error('Axon location generation failed to produce axons contained in fascicles'), 
    end
    n_ratio = n_ratio * 0.8; continue
  else xy = xy(1:nA,:); f_id = f_id(1:nA,:); vp = vp(1:nA); break
  end
end

uatt = upper(axon_types.Type);

isa_C_type = contains(uatt,'UNMY');
isa_aff = contains(uatt,'AFF') | contains(uatt,'SENS');

nAA = sum(axon_types.Count(~isa_C_type & isa_aff));
nAE = sum(axon_types.Count(~isa_C_type & ~isa_aff));
nCA = sum(axon_types.Count(isa_C_type & isa_aff));
nCE = sum(axon_types.Count(isa_C_type & ~isa_aff));

perimeter_fcn = @(v) sum(sqrt(sum((cp(v,:)-cp(v([2:end 1]),:)).^2,2)));
[~,is_myelin] = sort(cellfun(perimeter_fcn, vp),'descend');
is_myelin = ismember((1:nA)',is_myelin(1:nAA+nAE));

is_afferent = is_myelin;
[~,aidx] = sort(rand(nAA+nAE,1));
is_afferent(is_myelin) = ismember(1:sum(is_myelin),aidx(1:nAA));

[~,aidx] = sort(rand(nCA+nCE,1));
is_afferent(~is_myelin) = ismember(1:sum(~is_myelin),aidx(1:nCA));

if isfield(D,'fibreEfferent')
  error apply_type_customisation_Adelta
end

if isfield(D,'unmyelinatedEfferent')
  %%
  [~,aidx] = sort(rand(size(pop.unmyelinated_diam))); 
  pop.unmyelinated_diam = pop.unmyelinated_diam(aidx); 
  
  unmy_aff = is_afferent(~is_myelin);
  
  x = mean(diff(D.unmyelinatedEfferent.x))/2;
  x = D.unmyelinatedEfferent.x-x; 
  x(end+1) = inf; x(1) = -inf;
  
  for ii = 1:numel(x)-1
    
    sel = (pop.unmyelinated_diam >= x(ii) & ...
           pop.unmyelinated_diam < x(ii+1));
   if ~any(sel), continue, end
   unmy_aff(sel) = linspace(0,1,sum(sel)) < D.unmyelinatedEfferent.y(ii); 
  end
  
  del = sum(unmy_aff) - axon_types(1,2);
  
  if del > 0
    aidx = find(unmy_aff);
    [~,sel] = sort(rand(size(aidx))); 
    unmy_aff(aidx(sel(1:del))) = 0;
  elseif del < 0
    aidx = find(~unmy_aff);
    [~,sel] = sort(rand(size(aidx))); 
    unmy_aff(aidx(sel(1:del))) = 1;
  end
  
  is_afferent(~is_myelin) = unmy_aff; 
    
 %%
end

pop.fibre_xy = xy(is_myelin,:);
pop.fibre_fascicle = f_id(is_myelin);
pop.fibre_afferent = is_afferent(is_myelin);

pop.unmyelinated_xy = xy(~is_myelin,:);
pop.unmyelinated_fascicle = f_id(~is_myelin);
pop.unmyelinated_afferent = is_afferent(~is_myelin);

% TODO : code to duplicate fascicle1 configuration for replicated fascicles
% (I wrote this , but somehow lost it possibly in a git pull or an experiment folder)

% pop = patch_fascIDs(pop,nerve); 

assignin('caller','xy',xy);
assignin('caller','f_id',f_id);
assignin('caller','is_myelin',is_myelin);
assignin('caller','is_afferent',is_afferent);

return

%%

clf %#ok<UNRCH>
scatter(pop.fibre_xy(:,1),pop.fibre_xy(:,2),10*pop.fibre_diam,pop.fibre_fascicle,'o')
colormap(lines(max(pop.fibre_fascicle)))

%% are xy contained in the loop (alias inpolygon)
function is = in_loop(xy,loop), 
is = inpolygon(xy(:,1),xy(:,2),loop(:,1),loop(:,2));


%% Save each population as a .csv or .xlsx table    
function export_tables(P,sam,nerve,ext,dest)

if ~isfolder(dest), mkdir(dest); end
fprintf('Saving tables (.%s) to %s\n',ext,tools.file('T',dest))
for ty = 1:numel(P)
  %% Organise into tables 
  T = rmfield(P(ty),{'axon_type','axon_model','color'});
  T.myelinated = repmat(T.myelinated,size(T.fascicle));
  T.afferent = repmat(T.afferent,size(T.fascicle));
  
  T.axon_x = T.axon_xy(:,1);
  T.axon_y = T.axon_xy(:,2);
  T.axon_xy = []; 
  
  v = fieldnames(T); 
  v(~cellfun(@(x) isempty(T.(x)),v)) = [];
  T = struct2table( rmfield(T,v) );
  
  T = T(:,[1:2 end-(3:-1:0) 3:end-4]);
      
  output_file = tools.file('get','%s/axons %s (%s) (%%d).%s','next',dest,P(ty).axon_type,P(ty).axon_model,ext);
  writetable(T,output_file)

end

%%
T = struct;
T.ID = [1:numel(sam.A_diam) 1:numel(sam.C_diam)]'; 
T.myelinated = [true(size(sam.A_diam)) false(size(sam.C_diam))]';
T.diameter_um = [sam.A_diam sam.C_diam]';
T.g_ratio  = [sam.gratio ones(size(sam.C_diam))]';
T.axon_count = [arrayfun(@(x) sum(sam.A_axon == x),1:numel(sam.A_diam)) ...
                arrayfun(@(x) sum(sam.C_axon == x),1:numel(sam.C_diam))]';

T = struct2table(T);              

output_file = tools.file('get','%s/axon sampling (%%d).%s','next',dest,ext);  
writetable(T,output_file)

if isempty(nerve), return, end

nerve.axon_color = cat(1,P.color);

output_file = tools.file('get','%s/../anatomy/anatomy (%%d).%s','next',dest,'mat');
if ~isfolder(fileparts(output_file)), mkdir(fileparts(output_file)); end
fprintf('Saving %s\n',tools.file('T',output_file))
save(output_file, '-struct','nerve');

function nerve = get_anatomy(get_,named)


% Gather data or parse inputs for spatial layout and types 
if any(named('-anat'))    
  nerve = get_('-anat');
  if ischar(nerve)
    nerve = mesh.insert_gmsh_fascicles('-info','-file',nerve);
  else
    if ~isfield(nerve,'outline'), nerve = nerve.splines; end
    if ~isfield(nerve,'fascicles'), nerve.fascicles = nerve.outline; end
  end  
else nerve = mesh.insert_gmsh_fascicles('-info');  
end

if iscell(nerve.outline), nerve.fascicles = nerve.outline{1}; 
else nerve.fascicles = nerve.outline; 
end

if isfield(nerve,'Attributes') && isfield(nerve,'Children') % XML leftovers 
    nerve = rmfield(nerve,{'Name','Attributes','Data','Children','Type'});
end

if isfield(nerve,'splines'), nerve = rmfield(nerve,'splines'); end


