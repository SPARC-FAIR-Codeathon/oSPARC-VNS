
function [field,indices] = load_Ve_field(filename,varargin)
% This function supersedes load_S4L_field, and works for either
% EIDORS or SIM4LIFE extracellular potential (Ve) files. Not all EIDORS
% features are relevent if importing from S4L files. 
% 
% Options: 
%   -f [f_id] (default = fascicle 1)
%   -stim [1] : select stimulus ID from file. If -all is not set, generates
%               bipolar/tripolar etc stimuli assuming monopolar underlying
%               fields. 
%   -all : instead of compositing, return fields as a cell array with one
%          entry per -stim index. 
%   -pattern ... : set stimulus pattern, overrides -sim
%   -current ... : set current pattern matching -pattern. For
%                  instance, to specify a stimulus ganging together E1+2
%                  with a return of E3+4, use -pat 1:4 -current [1 1 -1 -1]
%                  defaults to [1 -1/n] = distributed return, appropriate
%                  for bipolar / tripolar / etc patterns.
% 
% Version 0.3 20-July-2020   Calvin Eiber refactor for BIDS model
% Version 0.2 16-April-2020  Calvin Eiber (defaults to "EMBC*" files)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};
contains = @(a,b) ~isempty(strfind(a,b)); %#ok<STREMP>

if nargin == 0
  filename = tools.file('get','sub~/eidors/stimulus*.mat','newest');
elseif isstruct(filename) % EIDORS data, pre-loaded
  
elseif strcmpi(filename,'S4L')
  filename = [expdir('embc*') 'HC_Sim4Life_fields.ve.dat']; 
elseif contains(filename,'~'), filename = tools.file(filename);
end


if ischar(filename) && strncmp(fliplr(filename),'tad.ev.',7) 
  %% ".ve.dat" = S4L file
  
  if any(named('-f')), fid = get_('-f'); % fascicle = ?
  else                 fid = 1;
  end
  
  dat = read_dat_file(filename); 
  axy = [dat.axons.x dat.axons.y]; 
  dat.centre = mean(axy);
  
  indices = (dat.axons.f == fid);
  [~,gx] = meshgrid(dat.axons.z,dat.axons.x(indices));
  [gz,gy] = meshgrid(dat.axons.z,dat.axons.y(indices));
  gz = gz*1000; % to units of mm 
  
  XYZ = [gx(:) gy(:) gz(:)];
  V = reshape(dat.voltage(indices,:),[],1);
  
else [XYZ,V,indices] = compose_from_eidors(filename); % EIDORS data
end


if ~isempty(which('scatteredInterpolant'))
  field = scatteredInterpolant(XYZ, V, 'natural');
  field.ExtrapolationMethod = 'nearest'; 
  if any(named('-all'))
    f = field; 
    field = cell(size(get_('-s'))); 
    field(:) = {f}; 
    for ii = 2:numel(field)
      [~,field{ii}.Values,~] = compose_from_eidors(filename,ii);
    end
  end
else
  % Define virtual v=0 corners:
  lim = [min(XYZ(:,3)) max(XYZ(:,3))];
  XYZ = [XYZ; lim([1 1 1; 1 1 2; 1 2 1; 2 1 1; 1 2 2; 2 1 2; 2 2 1; 2 2 2])];
  V = [V; 0;0;0;0; 0;0;0;0];
  % V(end+1 : size(XYZ,1)) = 0; 
  
  field = @(x,y,z) griddata3(XYZ(:,1),XYZ(:,2),XYZ(:,3),V,x,y,z); 
  
  % Wierd bug: if you get an error calling (e.g.) ffun(0,0,0)
  % error: griddatan: =: nonconformant arguments (op1 is 1x1, op2 is 4x1)
  % try calling with an extra dummy value in the x,y,z vector: seems to fix.
  
  if any(named('-all')) % inefficient but what can you do? 
    field = {field};     
    for ii = 2:numel(field)
      [~,V,~] = compose_from_eidors(filename,ii);
      field{ii} = @(x,y,z) griddata3(XYZ(:,1),XYZ(:,2),XYZ(:,3),V,x,y,z); 
    end
  end  
end


return

%% Compare results (not apples-to-apples)

if ~exist('vS4L','var'), XYZ_ = XYZ; vS4L = V; end

clf, subplot(2,1,1)
[~,pid] = max(vS4L);
pid = (XYZ_(:,3) == XYZ_(pid,3));
scatter(XYZ_(pid,1),XYZ_(pid,2),[],vS4L(pid),'o','filled')
axis equal, tools.tidyPlot

subplot(2,1,2), cla, hold on
[~,pid] = max((V));
pid = abs(XYZ(:,3) - XYZ(pid,3)) < 0.05
scatter3(XYZ(pid,1),XYZ(pid,2),XYZ(pid,3),[],V(pid),'o','filled')
axis equal, tools.tidyPlot




function [XYZ,V,indices] = compose_from_eidors(filename,stimIndex)

named = evalin('caller','named');
get_  = evalin('caller','get_');

if nargin < 2
 if any(named('-all')), stimIndex = 1; else stimIndex = []; end
end
     

if isstruct(filename), dat = filename;
else                   dat = load(filename);
end

if any(named('-f')), fid = get_('-f'); % fascicle = ?
else                 fid = 1;
end

current = []; 

if any(named('-pa')) % -pattern

  target = get_('-pa'); 
  stimuli = [dat.model.stimulation.stim_pattern];        
  if any(target == 0), current = target;
  else                 elec = target;
    nE = numel(elec); 
    if any(named('-c')), current = get_('-c'); 
    else current = [1 -ones(1,nE-1)/(nE-1)];
    end

    target = full(0*stimuli(:,1));
    target(elec) = current;      
    current = round(stimuli \ target,6);
  end    
  pid = 1:numel(current);
else % use what's present in the stimulus file

  if any(named('-s')), pid = get_('-s'); % sim # = ?
  else                 pid = 1;
  end
  if ~isempty(stimIndex), pid = pid(stimIndex); end
  if numel(pid) > 1 
    % If monopolar stimuli are what's in eidors\stimulus.mat, 
    % compositeto a multipolar stimulus. However, I think the default
    % stimulus is sequential bipoles:  E1.R2, E2.R3, etc.

    if any(named('-c')), current = get_('-c'); 
    else current = [1 -(1/(numel(pid)-1))];
    end
  end
end

indices = dat.model.object_id{strcmp(dat.model.object_name, ...
                                sprintf('Fascicle%d',fid))};
indices = unique(dat.model.elems(indices,:));
XYZ = dat.model.nodes(indices,[3 2 1]);
V   = dat.v_extracellular(indices,pid);  
if ~isempty(current), V = V * reshape(current,[],1); end

% cla, scatter3(XYZ(:,3),XYZ(:,2),XYZ(:,1),[],V), axis equal
