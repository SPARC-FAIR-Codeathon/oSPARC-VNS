
function preview_layout(e)
% Preview the layout of the electrode array defined in elec-geom.mat
% v0.1 Calvin Eiber 20-Apr-2020

if ~exist('e','var'), 
  if evalin('caller','exist(''e'',''var'')'), e = evalin('caller','e'); 
  else
    eidors_file = tools.parse_arguments({},'eidors','s*.mat'); 
    e = load(eidors_file,'info'); 
    e = e.info;
  end
end

nC = size(e.ElectrodePositions,1);
C = lines(nC);

clf, hold on, grid on

for ii = 1:nC
  
  xy = [-1 1 1 -1 -1; -1 -1 1 1 -1]/2 .* ...
       e.ElectrodeDimensions(e.ElectrodeTypeIndex(ii),[1 3])'; 
  plot(xy(1,:) + e.ElectrodePositions(ii,1), ...
       xy(2,:) + e.ElectrodePositions(ii,3), ...
       'Linewidth',1.8)
     text(e.ElectrodePositions(ii,1),e.ElectrodePositions(ii,3), ...
          num2str(ii),'Fontsize',12,'Color',C(ii,:),'Horiz','center')
end

axis equal tight, tools.tidyPlot
 
  
  


