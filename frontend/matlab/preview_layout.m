
function preview_layout(e,varargin)
% Preview the layout of the electrode array defined in elec-geom.mat
% options: -g greyscale (Default: lines colormap)
%          -t transpose xy
%          -o plot array outline (PDMS carrier)
% v0.1 Calvin Eiber 20-Apr-2020

if ~exist('e','var') || ~isstruct(e)
  if exist('e','var') && ~isstruct(e), varargin = [{e} varargin]; end
  if evalin('caller','exist(''e'',''var'')'), e = evalin('caller','e'); 
  else
    eidors_file = tools.parse_arguments({},'eidors','s*.mat'); 
    e = load(eidors_file,'info'); 
    e = e.info;
  end
end

named = @(v) strncmpi(v,varargin,length(v)); 

if isfield(e,'array'), e = e.array; end

nC = size(e.ElectrodePositions,1);
C = lines(nC);

if any(named('-g')), C(:) = 0.3; end

clf, hold on

if any(named('-t'))    
    e.ElectrodeDimensions = e.ElectrodeDimensions(:,[3 2 1]);
    e.ElectrodePositions = e.ElectrodePositions(:,[3 2 1]);
end

for ii = 1:nC
  
  xy = [-1 1 1 -1 -1; -1 -1 1 1 -1]/2 .* ...
       e.ElectrodeDimensions(e.ElectrodeTypeIndex(ii),[1 3])'; 
  plot(xy(1,:) + e.ElectrodePositions(ii,1), ...
       xy(2,:) + e.ElectrodePositions(ii,3), ...
       'Linewidth',1.8,'Color',C(ii,:))
     text(e.ElectrodePositions(ii,1),e.ElectrodePositions(ii,3), ...
          num2str(ii),'Fontsize',12,'Color',C(ii,:),'Horiz','center')
end


if any(named('-o')) % Draw outline
    
  G = @(v) [v v v]/10; 
  
  cir = [cos(linspace(0,2*pi,81)); sin(linspace(0,2*pi,81))]';
  if isfield(e,'CarrierOutline'), xy = e.CarrierOutline(:,[1 3]);
  else
      xy = [cir(1:21,:)*0.3  + [2.1 1.15]; ...
            cir(21:41,:)*0.3 - [2.1 -1.15]; 
            cir(41:61,:)*0.3 - [2.1 1.15]; ...
            cir(61:81,:)*0.3 + [2.1 -1.15]];
        xy = xy([1:end 1],:); 
  end
  
  if any(named('-t')), xy = xy(:,[2 1]); end
      
  plot(xy(:,1),xy(:,2),'k-','LineWidth',1.6,'Color',G(1.5)), hold on
  
  if isfield(e,'CarrierOutline2'), xy = e.CarrierOutline2(:,[1 3]);
  elseif ~isfield(e,'CarrierOutline')
      f_ = @flipud;
      xy = [ [ 1.5 2.9] + cir(1:41,:)*1; ...
            f_([0 2.9] - cir(1:41,:)/2); ...
             [-1.5 2.9] + cir(1:41,:)*1; ...            
            [-1.5 -2.9] + cir(41:81,:)*1; ...
            f_([0 -2.9] - cir(41:81,:)/2); ...
             [1.5 -2.9] + cir(41:81,:)*1];             
      xy( 42+12:82-12, 2) = 3.45; 
      xy = xy([1:end 1],:);   
  else xy = []; 
  end
  if ~isempty(xy)
    if any(named('-t')), xy = xy(:,[2 1]); end
    plot(xy(:,1),xy(:,2),'k-','LineWidth',1.0,'Color',G(1.5))
  end
end

axis equal tight, tools.tidyPlot, grid on

if isfield(e,'ElectodePatternID')
    title(e.ElectodePatternID,'Interpreter','none')
end


 
  
  


