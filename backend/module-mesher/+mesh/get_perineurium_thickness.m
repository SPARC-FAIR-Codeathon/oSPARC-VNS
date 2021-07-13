
function [widths, data] = get_perineurium_thickness(file)
% Get the perineurium thickness of each fascicle in the specified file

has_ext_ = @(a,b) strncmpi(fliplr(a),fliplr(b),length(b)); 

if nargin == 0, data = mesh.insert_gmsh_fascicles('-info'); 
    file = data.filename; 
elseif isstruct(file), data = file; 
elseif has_ext_(file,'.mat'), data = load(file); 
elseif has_ext_(source_file,'.xml') % from XML
  printf('Loading %s', source_file);
  source_anat = tools.parse_xml(source_file);
  source_anat.Type = 'MBF-XML';
  source_anat.filename = source_file;
else error('unknown filetype on file "%s", expected {.mat, .xml}', ...
                                     file)
end


if isfield(data,'outline') && iscell(data.outline)
  
    error todo_gather_correct_data
    widths = pw_from_outline(data.outline(1),data.outline(1))
    
elseif isfield(data,'Children') % most likely XML
    
    % XML handler microfunctions 
    the = @(x,n) x.Children(strncmpi({x.Children.Name},n,length(n)));
    attr_ = @(x,n) x.Attributes(strcmpi(n,{x.Attributes.Name})).Value;
    
    xy_ = @(x,n) arrayfun(@(p) str2double(p.Attributes(n).Value), x);
    
    profiles = the(data,'contour');
    names = arrayfun(@(p) attr_(p,'Name'),profiles,'unif',0);     
    fasc_idx = str2double(regexp(names,'(?<=ior-)\d+','match','once'));
    widths = zeros(max(fasc_idx),1);
    
    for ff = 1:max(fasc_idx)
        
        sel = (fasc_idx == ff); 
        
        xy_fasc = profiles(sel & contains(names,'Interior'));
        xy_peri = profiles(sel & contains(names,'Exterior'));
        
        if numel(xy_fasc) ~= 1 || numel(xy_peri) ~= 1
          warning('Fascicle #%d: found %d interior, %d exterior loops', ...
                             ff, numel(xy_fasc), numel(xy_peri))
          widths(ff) = NaN;
          continue
        end
        
        xy_fasc = the(xy_fasc,'point');
        xy_fasc = [xy_(xy_fasc,2); xy_(xy_fasc,3)]';
        
        xy_peri = the(xy_peri,'point');
        xy_peri = [xy_(xy_peri,2); xy_(xy_peri,3)]';
        
        widths(ff) = pw_from_outline(xy_fasc,xy_peri);
    end
else
    error('not sure how to parse the contents of "%s"', file)
end


return


function [mean_pw,pw] = pw_from_outline(fasc,peri)


DEBUG = false; 

% plot(peri.xy(:,1),peri.xy(:,2), ...
%      data(2).xy(:,1),data(2).xy(:,2), '-', 'LineWidth',1.2)
 
if isstruct(fasc), fasc = fasc.xy; end
if isstruct(peri), peri = peri.xy; end

pw = 0*peri(:,1); 
 
for ii = 1:size(fasc,1)    
    %%
    fxy = fasc(ii,:);    
    [~,idx] = min( (sum((peri - fxy).^2,2)) );     
    idx = mod(idx-[2 1 0],size(peri,1))+1;
    
    if DEBUG % Debug visualisations 
        clf, hold on
        plot(peri(idx,1),peri(idx,2),'-')
        plot(fxy(1),fxy(2),'o')
        axis equal
    end
       
    dhat = [1 -1 0; 0 -1 1] * peri(idx,:); 
    dhat = dhat ./ sqrt(sum(dhat.^2,2));
    
    uhat = (fxy-peri(idx(2),:)); 
    % uhat = uhat / sqrt(sum(uhat.^2)); 
    
    u = (dhat * uhat'); 
    d =  dhat.*max(0,u) + peri(idx(2),:);    
    
    pw(ii) = min(sqrt(sum((d - fxy).^2,2)));
    
    if DEBUG
        plot(d(:,1),d(:,2),'.') 
        plot(fxy(1)+pw(ii)*cos(linspace(0,2*pi,65)), ...
             fxy(2)+pw(ii)*sin(linspace(0,2*pi,65)),'-')
        % ginput(1); 
    end
    
    % plot(peri(idx(2),1) + [0 uhat(1)], peri(idx(2),2) + [0 uhat(2)], 'k--s')
    % plot(peri(idx(2),1) + dhat(1)*[0 u(1) 1], peri(idx(2),2) + dhat(3)*[0 u(1) 1], 'k-->')
    % plot(peri(idx(2),1) + dhat(2)*[0 u(2) 1], peri(idx(2),2) + dhat(4)*[0 u(2) 1], 'k--<')
    
    % d = [d(1,:); fxy; d(2,:)]; 
    % plot(d(:,1),d(:,2),':')    
    % pause(0.1)
end

mean_pw = mean(pw(2:end)); 