

if isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP>
  % Set up session defaults to match default MATLAB behaviour
  struct_levels_to_print(0)
  save_default_options('-mat-binary')
  more off
end