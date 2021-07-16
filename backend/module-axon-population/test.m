
if exist('compile','dir')
  disp('Testing COMPILED binaries')
  ! ./compile/for_testing/nerve_mesher rat-cervical-vagus "" ./input/demo/nerve-script.json 
  % ! ./compile/for_testing/nerve_mesher "./input/demo/array.json" "" "input/demo/nerve-script.json"
  return
end

%% EXAMPLE (1): planar-array

axons_file = 'rat-pelvic-nerve';
nerve_script = './input/demo/nerve-script.json';

make_axon_population(axons_file, ... % .mat file or selection from packaged files
              '',           ...  % as module_mesher
              nerve_script, ...  % as module_mesher
              'count', ... % one of {'count', 'per_mm2', 'ignore'}
              300,     ... % # myelinated axons 
              1000,    ... % # unmyelinated axons
              3.3)     % sensory : motor ratio ) 