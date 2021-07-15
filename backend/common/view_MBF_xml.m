
clearvars -except current_folder

if isempty(which('tools.parse_xml'))
    [tools_dir] = uigetdir('','Please locate the +tools folder to add to path for this session');    
    path(path, fileparts(tools_dir)); 
    fprintf('%s added to path\n', tools_dir); 
end

if ~exist('current_folder','var'), current_folder = '.'; end
[filename, current_folder] = uigetfile('*.xml', ... 
                                       'Select XML file to parse', ...
                                       current_folder);
fprintf('Loading %s from %s \n ... ', filename, current_folder)
xml = tools.parse_xml([current_folder filename]); 
fprintf('Done!\n')

% lambda functions to extract values
the = @(x,n) x.Children(strncmpi({x.Children.Name},n,length(n)));
attr_ = @(x,n) x.Attributes(strcmpi(n,{x.Attributes.Name})).Value;
trim_ = @(x) x.Children(~contains({x.Children.Name},'#text'));

%% Visualise and name counts 

contours = the(xml,'contour');

info.type = {}; 
info.count = []; 

clf, hold on

for ii = 1:numel(contours)
    
    this = contours(ii);     
    type = regexprep(attr_(this,'name'),'-\d+','');
    if ~ismember(type,info.type)
         info.type = [info.type {type}];
         info.count = [info.count 1]; 
    else sel = strcmpi(info.type,type); 
         info.count(sel) = info.count(sel)+1; 
    end

    points = the(this,'point'); 
    
    x = arrayfun(@(p) str2double(attr_(p,'x')), points);
    y = arrayfun(@(p) str2double(attr_(p,'y')), points); 
    
    try 
      c = attr_(this,'color');
      c = hex2dec([c(2:3);c(4:5);c(6:7)])'/255;        
    catch, c = [1 0 0];
    end
    
    plot(x,-y,'-','Color',c)
    text(x(1),-y(1),attr_(this,'name'),'FontSize',8,'Color',c)
    
end

axis equal, grid on, w = 0.03; 
axis(axis * [1+w -w 0 0; -w 1+w 0 0; 0 0 1+w -w; 0 0 -w 1+w]);

clear c w x y ii this type sel contours points 

arrayfun(@(a,b) fprintf('[%2d] %s\n', b, a{1}), info.type, info.count)




