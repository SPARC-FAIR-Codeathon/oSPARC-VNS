
if exist('compile','dir')
  disp('Testing COMPILED binaries')
  ! ./compile/for_testing/nerve_mesher rat-cervical-vagus "" ./input/demo/nerve-script.json per_mm2 2e4 7e4 3.3
  % ! ./compile/for_testing/nerve_mesher "./input/demo/array.json" "" "input/demo/nerve-script.json"
  return
end

%% EXAMPLE (1): "count"

axons_file = 'rat-pelvic-nerve';
nerve_script = './input/demo/nerve-script.json';

make_axon_population(axons_file, ... % .mat file or selection from packaged files
              '',           ...  % as module_mesher
              nerve_script, ...  % as module_mesher
              'per_mm2', ... % one of {'count', 'per_mm2', 'ignore'}
              2e4,     ... % # myelinated axons 
              7e4,    ... % # unmyelinated axons
              3.3)     % sensory : motor ratio ) 