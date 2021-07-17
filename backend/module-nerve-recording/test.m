
if exist('compile','dir')
    
    disp('Testing COMPILED binaries')
    
    ! ./compile/for_testing/module_EIDORS_compute_fields "./input/demo/example cuff.mat"    
    
    return
end
    
%% EXAMPLE (1): fields for cuff array

EIDORS_fwd_model; 

       
       