
if exist('compile','dir')
  disp('Testing COMPILED binaries')
  ! ./compile/for_redistribution_files_only/nerve_mesher ./input/C-FINE.json "" ./input/sub-57_sam-1.json
  
  % ! ./compile/for_testing/nerve_mesher "./input/demo/array.json" "" "input/demo/nerve-script.json"
  return
end

%% EXAMPLE (1): planar-array

tools.file('root',pwd)
nerve_mesh(tools.file('input~/demo/array.json'), '', ...
           tools.file('input~/demo/nerve-script.json'))

       
%% Preview set-up 2

opts = tools.parse_json('./input/demo/payne2019-nerve.json');
mesh.insert_gmsh_fascicles('-setup',opts)
plots.preview_fascicles

%% EXAMPLE (2): cuff-array

tools.file('root',pwd)
nerve_mesh(tools.file('input~/demo/payne2019-cuff.json'), '', ...
           tools.file('input~/demo/payne2019-nerve.json'))

       
%% EXAMPLE (3) - human vagus with epineurium

array = tools.parse_json( 'C-FINE.json' );
human_vagus = tools.parse_json( tools.file('get','sub*.json') );
mesh.insert_gmsh_fascicles('-setup',human_vagus)

% s = mesh.insert_gmsh_fascicles;
nerve = plots.preview_fascicles('-no-warn','-array',array.array);

%%

tools.file('root',pwd)
nerve_mesh('C-FINE.json', '', tools.file('get','sub*.json'))


           