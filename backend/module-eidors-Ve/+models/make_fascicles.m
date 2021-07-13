
function make_fascicles(nF, f_size_string, f_shape_string, f_roi_string, perineuri_string)
% models.make_fascicles : takes a series of input arguments describing the size
%    and shape distributions for a series of elliptical fascicles in a
%    rectangular or circular region of interest and use these 
% 
% # fascicles, fascicle_size, fascicle_shape, nerve_extent, perineurium
% 
% loop through size string to generate nF fascicles. If nF > number of
% elements in size string, loop back to the beginning. 
% 
%  elements are seperated with ';'
%    mean +- sd   : normally-distributed
%    LB-UB        : uniformly-distributed
% 
%   N~(mu,sigma)  : normally distributed (see makedist for definitions)
%   LN~(mu,sigma) : log-normally distributed 
%   L~(mu,sigma)  : logistically distributed
%   LL~(mu,sigma) : log-logistically distributed
%   E~(lambda)    : exponentially distributed
%   U~(LB, UB)    : uniformly distributed
%   W~(a,b)       : Weibull distribution

if nargin < 5, perineuri_string = '3%'; end % standard approximation
if nargin < 4, f_roi_string = '1'; end % in mm
if nargin < 3, f_shape_string = '0.6-0'; end % 
if nargin < 2, f_size_string = 'LL~(5.2,0.285)'; end % human vagus
if nargin < 1, nF = round(rand*10)+5; end
if ischar(nF), nF = str2double(nF); end

% human vagus based on dataset: Pelot, NA, et al (2020). 
% Quantified Morphology of the Human Vagus Nerve with Anti-Claudin-1 
% https://doi.org/10.26275/NLUU-1EWS

% the Pelot Vagus perineurium equation is 3.702%+10.5

outline = parse_outline_parameters(nF, f_size_string, ...
                                       f_shape_string, ...
                                       f_roi_string, ...
                                       perineuri_string);

%%

fi = fopen('output.xml','wt'); 
w_ = @(s,varargin) fprintf(fi,[s '\n'],varargin{:}); 

w_('<?xml version="1.0" encoding="ISO-8859-1"?>');
w_('<mbf version="4.0" xmlns="http://www.mbfbioscience.com" xmlns:nl="http://www.mbfbioscience.com" appname="Neurolucida 360" appversion="2020.1.1 (64-bit)">');
w_('<description><![CDATA[]]></description>');

% w_('<images>');
% w_('  <image>');
% w_('    <filename>%s</filename>',[pwd 'output.xml']);
% w_('    <scale x="%0.6f" y="%0.6f"/>',1,1);
% w_('  </image>');
% w_('</images>');

w_('<parameters>');
w_('  <parameter name="number of fascicles" value="%g"/>', nF);
w_('  <parameter name="nerve region" value="%s"/>', f_roi_string);
w_('  <parameter name="fascicle sizes" value="%s"/>', f_size_string);
w_('  <parameter name="fascicle eccentricity" value="%s"/>', f_shape_string);
w_('  <parameter name="perineurium thickness" value="%s"/>', perineuri_string);
w_('</parameters>');

type_codes = unique([outline.type]); 
type_index = ones(size(type_codes));

for oid = 1:numel(outline)
    
    type_id = (type_codes == outline(oid).type);
    
    if outline(oid).type == 'f', nom = 'FascicleInterior';
        meta = 'color="#4DBEEE" closed="true" shape="Contour"';
    elseif outline(oid).type == 'p', nom = 'FascicleExterior';
        meta = 'color="#EDB120" closed="true" shape="Contour"';
    elseif outline(oid).type == 'x', nom = 'NerveOutline';
        meta = 'color="#CCCCCC" closed="true" shape="Contour"';
    else error object_class_name
    end
    
    w_('<contour name="%s-%d" %s>', nom, type_index(type_id), meta) ;
    w_('  <units>microm</units>');
    
    for pp = 1:size(outline(oid).xy,1)
      w_('  <point x="%0.4f" y="%0.4f" z="0.00" d="0.01"/>', ...
                    outline(oid).xy(pp,:) .* [1 -1]);
    end      
    w_('</contour>');
end

w_('</mbf>');

fclose(fi); 
%%


function obj = parse_outline_parameters(nF, sSizes,sEccen,sROI,sPeri)


roi = str2double(regexp(sROI,'(\d*\.)?\d+','match')); 

switch(numel(roi))
  case 1, roi = [0 0 roi];
  case 2, roi = [-roi(1)/2 0 roi];
end

theta = linspace(0,2*pi,61)';

switch(numel(roi))
    
  case 3, bounds = roi(3)/2 * [cos(theta) sin(theta)];
          bounds = bounds + roi(1:2);
  case 4, bounds = [1 1 0 0 1; 0 1 1 0 0]'.*roi(3:4) + roi(1:2);
  otherwise    
    error('Unknown number of elements in "%s" (found %d, expected 1-4)', ...
                                        sROI, numel(roi))
end

bounds = 1e3 * bounds; % convert bounds to um

%% Extract distribtuion specifications for fascicle sizes 

sSizes = strsplit(sSizes,';');
sEccen = strsplit(sEccen,';');

size_fcn = [sSizes sEccen]; 

for ii = 1:numel(size_fcn)
    
    this = struct;
    this.str = size_fcn{ii};
    this.par = str2double(regexp(this.str,'(\d*\.)?\d+','match'));
    
    if any(ismember(this.str,'~()')) % dist~(pars)
      this.dist = regexp(this.str,'^[^~\(]*','match','once');
    elseif any(ismember(this.str,'+�'))
      this.dist = 'N'; 
    else this.dist = 'U'; 
    end
    
    this.dist(ismember(this.dist,'a':'z')) = ''; % remove lowercase    
    if strcmp(this.dist,'U'), this.par = sort(this.par); end
        
    switch this.dist
        case 'E', pd = makedist('Exponential','mu',this.par(1));
        case 'L', pd = makedist('Logistic','mu',this.par(1),'sigma',this.par(2));
        case 'N', pd = makedist('Normal','mu',this.par(1),'sigma',this.par(2));
        case 'U', pd = makedist('Uniform','lower',this.par(1),'upper',this.par(2));
        case 'W', pd = makedist('Weibull','a',this.par(1),'b',this.par(2));
        case 'LN', pd = makedist('LogNormal','mu',this.par(1),'sigma',this.par(2));
        case 'LL', pd = makedist('LogLogistic','mu',this.par(1),'sigma',this.par(2));
        otherwise error('unknown distribution "%s", expected one of {%s}', ...
                         this.dist, 'N,LN,L,LL,U,E,W')
    end

    this.dist = pd;
    size_fcn{ii} = this;
end

eccen_fcn = size_fcn( (end-numel(sEccen)+1) : end);
size_fcn  = size_fcn( 1:numel(sSizes) );

s_peri = strrep(sPeri,'%','/100*x');
peri_fun = str2func(['@(x) ' s_peri]); 

%%

ff = 0; 
% v_ = @(x) reshape(x,[],6);
xy_ = @(t,ab,y) ab.*[cos(t+y) sin(t+y)]*[cos(y) -sin(y); sin(y) cos(y)];

fail_counter = 0; 
obj = [];

inpoly = @(q,b) inpolygon(q(:,1),q(:,2),b(:,1),b(:,2));

display_steps = false;


while ff < nF
    
    
    ff = ff + 1; 
    
    if fail_counter > 999
      error('Fascicle packing aborted at fascicle %d after %d failed iterations', ff, fail_counter)
    end
    
    
    d = mod(ff-1,numel(size_fcn)) + 1;
    e = mod(ff-1, numel(eccen_fcn)) + 1; 
    
    %% Get diameter / eccentricity / angle from CDFs 
    
    diam = size_fcn{d}.dist.icdf(rand);
    ecc  = eccen_fcn{e}.dist.icdf(rand);
    dpp = diam + peri_fun(diam); % diameter + perineurium
    ang  = 2*pi*rand;         
    
    a = dpp/2 * (1-ecc.^2).^(0.25);
    b = dpp/2 / (1-ecc.^2).^(0.25);
    
    xy_outline = xy_(theta,[a b],ang);
  
%%
    
    if numel(bounds) == 10
         xy_mm = max(xy_outline);
         xy_b = bounds - (xy_mm .* [ 1 -1; 1 1; -1 1; -1 -1; 1 -1]);
         xy_b(end,:) = []; 
    else xy_b = bounds - xy_outline;        
    end
    
    for ii = 1:2:numel(obj)
      % plot(obj(ii).xy(:,1),obj(ii).xy(:,2),'b-')
      % plot(obj(ii).xy(1,1),obj(ii).xy(1,2),'b.')
      % plot(obj(ii).xy(:,1)+xy_outline(:,1), ...
      %      obj(ii).xy(:,2)+xy_outline(:,2),'g.')
      
      xy_b = [xy_b; obj(ii).xy(2:end,:) + xy_outline(2:end,:)];
    end
    
    ok = inpoly(xy_b,bounds);
%     for ii = 1:2:numel(obj)
%       ok = ok & ~inpoly(xy_b,obj(ii).xy + xy_outline);
%     end
    
    if ~any(ok)
       ff = ff - 1; 
       fail_counter = fail_counter + 1;
       continue
    end
    
    xy_b = xy_b(ok,:); 
       
    %%
    
    tri = delaunay(xy_b);
    cxy = squeeze(mean(reshape(xy_b(tri,:),[],3,2),2));
    ok = inpoly(cxy,bounds);
        
    for ii = 1:2:numel(obj)
        ok = ok & ~inpoly(cxy,obj(ii).xy + xy_outline);
    end
    
    tri = tri(ok,:); cxy = cxy(ok,:); %#ok<NASGU>
    % plot(cxy(ok,1),cxy(ok,2),'.','Color',[.7 .7 .7])
    
    if ~any(ok)
       ff = ff - 1; 
       fail_counter = fail_counter + 1;
       continue
    end 
    
    abc = sqrt( (xy_b(tri,1) - xy_b(tri(:,[2 3 1]),1)).^2 + ...
                (xy_b(tri,2) - xy_b(tri(:,[2 3 1]),2)).^2 );
    abc = reshape(abc,[],3);
    
    % where s is the semi-perimeter of the triangle
    s = sum(abc,2)/2;
    t_a = real(sqrt(s.*prod(s-abc,2))); % area (Herod's formula)
    
    % validation on rand() approach used here: 
    % PDF = [2 4 0 3 1];
    % for ii = 1:1e5, sel(ii) = find(rand*sum(PDF) < cumsum(PDF),1); end    
    
    sel = find(rand*sum(t_a) < cumsum(t_a),1); % area = p(pick)
    [~,seq] = sort(rand(1,3));  
    txy = xy_b(tri(sel,seq),:); % shuffle order of vertices 
    u = rand(2,1); u(1) = u(1).^2;
    uxy = u(1)*(txy(2,:)-txy(1,:)) + ...
          u(2)*(1-u(1))*(txy(3,:)-txy(1,:)) + txy(1,:);
    
    if display_steps
        clf, hold on, axis equal %%#ok<UNRCH>
        plot(bounds(:,1),bounds(:,2),'k-',xy_b(:,1),xy_b(:,2),'r.')
        plot(cxy(:,1),cxy(:,2),'.','Color',[.6 .6 .6])

        plot(txy([1:end 1],1),txy([1:end 1],2),'r-')
        plot(uxy(1),uxy(2),'b+')
        plot(xy_outline(:,1)+uxy(1), xy_outline(:,2)+uxy(2),'b-')
        % ginput(1);
    end
    
    %%
    
    a = diam/2 * (1-ecc.^2).^(0.25);
    b = diam/2 / (1-ecc.^2).^(0.25);
        
    this = struct;
    this.type = 'p';
    this.xy = xy_outline + uxy;
    
    if ~all(inpoly(this.xy,bounds))
         ff = ff - 1; 
         fail_counter = fail_counter + 1;
         continue
    else fail_counter = 0; 
    end
    
    
    if isempty(obj), obj = this;
    else obj(end+1) = this;
    end
    
    this = struct;
    this.type = 'f';
    this.xy = xy_(theta,[a b],ang) + uxy;
    obj(end+1) = this;
    
    %%
end

if 1
    %% Post-generation display
    clf, hold on, axis equal %%#ok<UNRCH>
    plot(bounds(:,1),bounds(:,2),'-','Color',[0 0 0 0.3])    
    C = lines(numel(obj));
    for ii = 1:numel(obj)
      if mod(ii,2)
           style = {'-','Color',(C((ii+1)/2,:)+1.2)/2.2,'LineWidth',1.4};
      else style = {'-','Color',C(ii/2,:),'LineWidth',1.0};
          
          cxy = mean(obj(ii).xy(2:end,:));
          text(cxy(1),cxy(2),num2str(ii/2),style{2:3},'FontSize',7)
      end
      plot(obj(ii).xy(:,1), obj(ii).xy(:,2),style{:})
      
    end
        
    tx = unique( [1;-1]*(0:100:max(abs([xlim ylim]))) );     
    w = 0.02;
    
    tools.tidyPlot, axis tight, grid on
    axis(axis*[1+w -w 0 0; -w 1+w 0 0; 0 0 1+w -w; 0 0 -w 1+w])
    set(gca,'XTick',tx,'YTick',tx,'XColor','none','YColor','none')
    
    plot(tx(end-[1 2]), (ylim*[1.02;-0.02])*[1 1], 'k-','LineWidth',1.3,'Clipping','off')
    text(mean(tx(end-[1 2])), ylim*[1.03;-0.03], '100 \mum', 'Horiz','center','vert','top')    
end

%%


this.type = 'x';
this.xy = bounds; 
obj(end+1) = this;


return



