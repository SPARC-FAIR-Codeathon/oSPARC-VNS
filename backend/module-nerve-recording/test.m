
if exist('compile','dir')
    
    disp('Testing COMPILED binaries')
    
    system(sprintf('%s "%s" "%s" "%s" "" 24.411', ...
                   './compile/for_redistribution_files_only/simulate_nerve_recording', ...
                   './input/demo/axon-population (1).mat', ...
                   './input/demo/extracellular-potential (1).mat', ...
                   'burst[0.2,20,5]_r2'));
    
    % ! ./compile/for_testing/simulate_nerve_recording "" "" "" "" 30
    
    return
end
    
%% EXAMPLE (1): fields for cuff array

make_nerve_recording('./input/demo-2/axon-population.mat', ...
                     './input/demo-2/extracellular-potential (1).mat','flat[1]_r1','')