
if exist('compile','dir')
  disp('Testing COMPILED binaries')
  ! ./compile/for_testing/nerve_mesher "./input/demo/array.json" "" "input/demo/nerve-script.json"
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

       
%% EXAMPLE (3)

human_vagus = tools.parse_json( tools.file('get','*.json') );
mesh.insert_gmsh_fascicles('-setup',human_vagus)
nerve = plots.preview_fascicles;

array = tools.parse_json('./input/demo/payne2019-cuff.json');

hold on, plot([1 -1 -1 1 1]*array.array.carrier.cuff_IDx/0.002, ...
              [1 1 -1 -1 1]*array.array.carrier.cuff_IDy/0.002, 'k-', ...
              'LineWidth', 1.5)
          

%%

tools.file('root',pwd)
nerve_mesh(tools.file('input~/demo/payne2019-cuff.json'), ...
                        human_vagus_sample)


           