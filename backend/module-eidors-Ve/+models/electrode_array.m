function e = electrode_array(e,pattern_ID)
% Generate electrode array
% 
% Typically I use this to generate a pre-specified array, but this can also
% parse an input structure to generate an electrode array pattern spec
% which can be digested by mesh.insert_gmsh_electrodes. 
% 
% The following fields will be generated if not present at input: 
%   e.DomainSize : ±[x y z]
%   e.ElectrodePositions  : [x y z], where y is normal to the planar array 
%                                      and the nerve lies along the x-axis
%   e.ElectrodeDimensions : [x (width) y_inset z (height)]
%   e.ElectrodeTypeIndex  : for multiple rows of e.ElectrodeDimensions,
%                            what size/shape is the nth electrode?
% 
% the following fields will be digested if present:
%   BipolarWithinPairSpacing
%   BipolarBetweenPairSpacing, BetweenElectrodeSpacing
%   n_ElectrodePairs (generates bipolar electrodes)
%   n_Tripoles   
%   n_Electrodes (generates equispaced monopolar electrodes
%   ElectrodeSize (alias for ElectrodeDimensions)
% 
% if a non-standard array is generated, plots.preview_layout will be called
%   to make sure what you've generated is what you wanted. 


if nargin == 1 
  if ~isstruct(e), pattern_ID = e; e = struct; 
  elseif isfield(e,'ElectodePatternID'), pattern_ID = e.ElectodePatternID;
  elseif isfield(e,'Name'), pattern_ID = e.Name;
  else pattern_ID = inputname(1);
  end
elseif nargin == 0, pattern_ID = 'A'; e = struct; 
end

switch pattern_ID
  %% Bionics Institute patterns
  case 'A'
    e.ElectrodePositions = [-1.85 0 0; -1.10 0 0; ...
                             1.10 0 0;  1.85 0 0];
    e.ElectrodeDimensions = [ 0.2 0.1 1.8 ];
    e.ElectrodeTypeIndex  = [ 1 1 1 1 ]; 
    e.ElectodePatternID = 'Bionics Institute "SPARC A"';

  case 'B'
    e.ElectrodePositions = [ -2 0 0; -1.25 0 0; ...
                           1.25 0 0;  2    0 0];
    e.ElectrodeDimensions = [ 0.2 0.1 1.8 ];
    e.ElectrodeTypeIndex  = [ 1 1 1 1 ]; 
    e.ElectodePatternID = 'Bionics Institute "SPARC B"';
    
    otherwise    
    %% Generate electrode layout pattern
    v_ = @(x) reshape(x,[],1);
    
    if isfield(e,'ElectrodeSize'), e.ElectrodeDimensions = e.ElectrodeSize;
        % e = rmfield(e,'ElectrodeSize');
    elseif ~isfield(e,'ElectrodeDimensions')
        e.ElectrodeDimensions = [0.2 1.8];
    end
    
    if size(e.ElectrodeDimensions,2) == 2
        e.ElectrodeDimensions = e.ElectrodeDimensions(:,[1 1 2]) ...
                                                .* [1 0 1] + [0 0.1 0];
    end
    
    if isfield(e,'BipolarWithinPairSpacing'), 
            w = e.BipolarWithinPairSpacing;
    else    w = 0.75; 
    end
    
    if isfield(e,'BipolarBetweenPairSpacing'), 
         b = e.BipolarBetweenPairSpacing;
    elseif isfield(e,'BetweenElectrodeSpacing'), 
         b = e.BetweenElectrodeSpacing;
    else b = 3.25;     
    end
    
    if ~isfield(e,'ElectrodePositions')        
      if isfield(e,'n_ElectrodePairs')
        e.ElectrodePositions = (1:e.n_ElectrodePairs)' * [b 0 0];
        e.ElectrodePositions = [e.ElectrodePositions + [w/2 0 0]; ...
                                e.ElectrodePositions - [w/2 0 0]];
        e.ElectrodePositions(:,1) = sort(e.ElectrodePositions(:,1)) - ...
                                    mean(e.ElectrodePositions(:,1));
      elseif isfield(e,'n_Bipoles')
        e.ElectrodePositions = (1:e.n_Bipoles)' * [b 0 0];
        e.ElectrodePositions = [e.ElectrodePositions + [w/2 0 0]; ...
                                e.ElectrodePositions - [w/2 0 0]];
        e.ElectrodePositions(:,1) = sort(e.ElectrodePositions(:,1)) - ...
                                    mean(e.ElectrodePositions(:,1));
      elseif isfield(e,'n_Electrodes')
        e.ElectrodePositions = (1:e.n_Electrodes)' * [b 0 0];
        e.ElectrodePositions(:,1) = sort(e.ElectrodePositions(:,1)) - ...
                                    mean(e.ElectrodePositions(:,1));  
      elseif isfield(e,'n_Tripoles')
        e.ElectrodePositions = (1:e.n_Tripoles)' * [b 0 0];
        
        e.ElectrodePositions = [e.ElectrodePositions + [w 0 0]; ...
                                e.ElectrodePositions; ...
                                e.ElectrodePositions - [w 0 0]];
        e.ElectrodePositions(:,1) = sort(e.ElectrodePositions(:,1)) - ...
                                    mean(e.ElectrodePositions(:,1));        
      else e.ElectrodePositions = v_([-b-w -b+w  b-w b+w]/2) * [1 0 0];
      end
    elseif size(e.ElectrodePositions,2) == 1
      e.ElectrodePositions = e.ElectrodePositions * [1 0 0];
    elseif size(e.ElectrodePositions,2) == 2
      e.ElectrodePositions = e.ElectrodePositions(:,[1 2 2]) .* [1 0 1];
    elseif size(e.ElectrodePositions,2) > 3
        e.ElectrodePositions = v_(e.ElectrodePositions) * [1 0 0];
    end
    
    
    e.ElectrodeTypeIndex = mod(0:size(e.ElectrodePositions,1)-1, ...
                                 size(e.ElectrodeDimensions,1)) + 1;
    e.ElectodePatternID = pattern_ID;
    
    plots.preview_layout(e)
end

e.ElectodePatternID = pattern_ID;

