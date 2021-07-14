
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
