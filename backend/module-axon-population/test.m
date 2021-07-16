
if exist('compile','dir')
  disp('Testing COMPILED binaries')
  ! ./compile/for_testing/generate_axon_population rat-cervical-vagus "" ./input/demo/nerve-script.json per_mm2 4e4 1.5e5 3.3  
  return
end

%% EXAMPLE (1): "count"

axons_file = 'rat-pelvic-nerve';
nerve_script = './input/demo/nerve-script.json';

make_axon_population(axons_file, ... % .mat file or selection from packaged files
              '',           ...  % as module_mesher
              nerve_script, ...  % as module_mesher
              'per_mm2', ... % one of {'count', 'per_mm2', 'ignore'}
              4e4,     ... % # myelinated axons 
              1.5e5,    ... % # unmyelinated axons
              3.3)     % sensory : motor ratio ) 