
if exist('compile','dir')
  disp('Testing COMPILED binaries')
  s_ = @(varargin) system(sprintf(varargin{:}));
  
  s_('%s "%s" "%s" "%s" per_mm2 40000 150000 3.3', ...
     './compile/for_testing/generate_axon_population', ...
     './input/rat-cervical-vagus.mat', ...
     '../module-mesher/input/sub-57_sam-1.xml', ...
     '../module-mesher/input/sub-57_sam-1.json')
  
  % ! ./compile/for_testing/generate_axon_population rat-cervical-vagus "" ./input/demo/nerve-script.json per_mm2 4e4 1.5e5 3.3  
  return
end

%% EXAMPLE (1): "per_mm2"

make_axon_population('./input/rat-cervical-vagus.mat', ... 
                     '../module-mesher/input/sub-57_sam-1.xml', ...  % as module_mesher
                     '../module-mesher/input/sub-57_sam-1.json', ...  % as module_mesher
                     'per_mm2', ... % one of {'count', 'per_mm2', 'ignore'}
                      4e4,     ... % # myelinated axons 
                      1.5e5,    ... % # unmyelinated axons
                      3.3)     % sensory : motor ratio ) 