

function makeFromTemplate(template_file,varargin)
% makeFromTemplate(template_file, ... )
% 
% makeFromTemplate reads a *.*.template file and evaluates meta-
%    expressions in that file, outputting a file with name "*.*". This
%    allows Matlab to modify model scripts for various other pieces of
%    software (GMSH, Neuron) while enabling the underlying scripts to be
%    written mostly in the target interpreter's source language with only a
%    few meta-expressions to allow matlab to, e.g. perform parameter
%    sweeps. 
% 
% Currently, three kinds of meta-expression are implemented: 
%   #ifdef / #endif blocks
%   #postpro (#post, #postprocess) blocks
%   $variable-name replacement
%   @[expression] evaluation
% 
% if the .*.template file contains #ifdef X / #endif, then anything which 
%   falls between those lines is only included if the string "X" was passed
%   into makeFromTemplate. 
% 
% if the .*.template file contains $variable-name expressions, the string
%   "$X" is replaced with the argument "X" (passed in using key-value
%   syntax e.g. make_( template , 'radius', 5) OR with the 
%   matlab variable "X" in the local workspace. The expression 
%    >> make_( template, 'radius', 5); is equivalent to
%    >> radius = 5; make_( template );
% 
% generate_gmsh_mesh will generate a helpful error if a $variable is used 
%   in the .*.template file which can't be found. A double "$$" is used to
%   indicate a literal "$" in the output file; no variable replacement is
%   performed. 
% 
% You can also call >> make_(file,'-makev')) which will create a file with
%   extension .*.variable, which can then be used to define variables and
%   parameter sweeps using the syntax >> make_(file,'-usev'); or
%    >> make_(file,'-sweep'); for parameter sweeps. 
% 
% if the .*.template file contains @[expression], the matlab code 
%   "expression" will be evaluated in the local workspace. This can be
%   combined with variable names 
%     (e.g. "Mesh.CharacteristicLengthMax = @[$radius / 100];" )
%   or used to evaluate arbitrary matlab code 
%     (e.g. "Mesh.CharacteristicLengthMax = @[my_lenmax_fun(args)];" )
% 
% For both $var and @[expr] replacement, numerical values are converted to
%   text, and text is used as-is. Attempting to insert non-numeric, non-
%   textual information will result in a lightly corrupted output file. 
% 
% By default, the output file has the same name as the input file but a
%   different output file can be specified e.g.
% 
%   >> _make(file,'-output',output_filename)
% 
% By: Calvin Eiber (calvin.eiber@unimelb.edu.au)


named = @(v) strncmpi(v,varargin,length(v)); 

if nargin == 0 || isempty(template_file), % prompt the user for a file    
   [template_file,t_path] = uigetfile('*.template');
    template_file = [t_path template_file];
end

if any(named('-output'))
     output_file = varargin{find(named('-output'))+1};
else output_file = strrep(template_file,'.template','');
end

if any(named('-makev')) % make variable file from template
    make_varfile_from_template(template_file, varargin{:})
    return
end

input_key = {}; 
input_val = {}; 
input_def = {}; 
define_key = {};
define_val = []; 

if any(named('-usev')) || any(named('-sweep'))
  [input_key,  input_val, input_def] = import_from_variable_file(template_file,varargin{:});
  for ii = 1:length(input_key)
    if ~any(named(input_key{ii})), continue, end
    input_val(ii) = varargin(find(named(input_key{ii}))+1); 
  end
  if any(named('-sweep')), return, end % escape recursion 
end
  
%% Make file from template 

t = fopen(template_file,'rt'); 
t_closer = onCleanup(@() fclose(t)); 

if t <= 0 % Multiple processes trying to read the same resource
 for n_try = 1:20
   if isempty(getCurrentTask), break, end
   pause(rand/20)
   t = fopen(template_file,'rt');
   if t > 0, break, end
 end
 assert(t > 0, 'failed to open template file');
end

g = fopen(output_file,'wt');
g_closer = onCleanup(@() fclose(g)); 

postprocess = {}; 

while ~feof(t)
  
  in = fgetl(t); 
  if ~ischar(in), continue, end
  
  postprocess_this_line = false; 
  
  if isempty(in), % Blank line
    if all(define_val), fprintf(g,'\n'); end
    continue, 
  end
  
  while strncmp(fliplr(strtrim(in)),'...',3)
    in = [regexprep(in,'\.{3}.*','') ' ' strtrim(fgetl(t))];
  end
  
  if strncmp(strtrim(in),'#IFDEF',6)
    define_key{end+1} = regexp(in,'(?<=#IFDEF[\s]+)\w+','match','once'); %#ok<AGROW>
    define_val(end+1) = any(named(define_key{end})) || ...
                        any(strcmp(define_key{end}, input_def)); %#ok<AGROW>
    fprintf(g,'// %s = %d\n',in,define_val(end));
    continue
  end
  
  if strncmp(strtrim(in),'#ENDIF',6)
    assert(~isempty(define_key),'Improperly formatted template')
    define_key(end) = []; 
    define_val(end) = []; 
    continue
  end
  
  if ~all(define_val), continue, end 
  
  % debug_check = any(ismember(in,'[]'));
  % if debug_check, disp(in), end
  
  % parse $variable meta-expressions
  tokens = unique(regexp(in,'(?<=[^\$]\$)\w+','match')); 
  for v = 1:numel(tokens), tok = tokens{v}; 
    
    if any(strcmp(input_key,tok)), k = strcmp(input_key,tok);
      if isnumeric(input_val{k})
        input_val{k} = num2str(input_val{k},12); %#ok<AGROW>
      end
      in = regexprep(in,['(?<=[^\$])\$' tok '([^\w])'],[input_val{k} '$1']);
      continue
    elseif any(named(tok))
      input_key{end+1} = tok; %#ok<AGROW>
      input_val{end+1} = varargin{find(named(tok))+1}; %#ok<AGROW>
    elseif evalin('caller',['exist(''' tok ''',''var'')'])      
      input_key{end+1} = tok; %#ok<AGROW>
      input_val{end+1} = evalin('caller',tok); %#ok<AGROW>
    else error('Undefined token "%s" in template file', tok)
    end
    
    if isnumeric(input_val{end})
      input_val{end} = num2str(input_val{end},12); 
    end
    in = regexprep(in,['\$' tok '([^\w]|$)'],[input_val{end} '$1']);
  end
    
  
  if strncmp(strtrim(in),'#POST',5)
    postprocess{end+1} = regexprep(strtrim(in),'#POST[^\s]*\s*','');  %#ok<AGROW>
    continue
  end
  
  
  % parse @[ ... ] meta-expressions
  tokens = unique(regexp(in,'\@\[[^\]]+\]','match')); 
  for v = 1:numel(tokens) % evalin(caller) 
    expr = tokens{v}(3:end-1);
    expr = evalin('caller',expr);
    if isnumeric(expr)
      expr = num2str(expr,12); 
    end % convert if needed
    in = strrep(in,tokens{v},expr);
  end
  in  = strrep(in,'$$','$'); % handle escaped "$" 
  
  fprintf(g,'%s\n',in);
  
  % if debug_check, disp(in), end
end

% fclose(t);
% fclose(g); 
clear t_closer g_closer 
for pp = 1:numel(postprocess)
  evalin('caller',[postprocess{pp} ';']);
end


  
  
  
%%


function make_varfile_from_template(template_file,varargin)

output_file = strrep(template_file,'.template','.variable');
% named = @(v) strncmpi(v,varargin,length(v)); 

%% Read defined #IFDEF and $variable strings from .template file

t = fopen(template_file,'rt'); 
assert(t > 0, 'failed to open template file');
t_closer = onCleanup(@() fclose(t)); 

def_name = {};  
def_context = {}; 
var_name = {}; 
var_context = {}; 

line_no = 0; 

while ~feof(t)
  
  in = fgetl(t);  line_no = line_no+1;
  if ~ischar(in), continue, end
  if isempty(in), continue, end % Blank line 
  
  while ~isempty(regexp(strtrim(in),'\.{3}','once')) % "..."
    in = [regexprep(in,'\.{3}.*','') ' ' strtrim(fgetl(t))];
    line_no = line_no+1;
  end
  
  if strncmp(strtrim(in),'#IFDEF',6)      
    def_name{end+1} = regexp(in,'(?<=#IFDEF[\s]+)\w+','match','once'); %#ok<AGROW>    
    def_context{end+1} = sprintf('// #%d: %s', line_no,in); %#ok<AGROW>
    continue
  end
  
  % parse $variable meta-expressions
  tokens = unique(regexp(in,'(?<=[^\$]\$)\w+','match')); 
  for v = 1:numel(tokens)
      var_name(end+1) = tokens(v);                            %#ok<AGROW>
      var_context{end+1} = sprintf('// #%d: %s', line_no,in); %#ok<AGROW>
  end
end

clear t_closer % fclose(t); 

%%

v = fopen(output_file,'wt');
v_closer = onCleanup(@() fclose(v)); 
assert(v > 0, 'failed to open variable file');
fprintf(v,'\n// Variable definition file generated using %s\n//  at %s\n', mfilename, datestr(now));

% Do #define statements first
[def_name,first_,indices] = unique(def_name,'stable');
fprintf(v,'\n#DEFINE statements\n');

for ii = 1:length(def_name)
    
    context = def_context{first_(ii)};
    
    if sum(indices == ii) > 1        
      if length(context) > 55, 
           context = sprintf('%s ... (used %d times)',context(1:50), ...
                                                        sum(indices == ii));
      else context = sprintf('%s (used %d times)',context, ...
                                                       sum(indices == ii));
      end
    elseif length(context) > 75, context = [context(1:70) ' ... ']; 
    end
    
    fprintf(v,'\n%s\n%s\n', context, def_name{ii});
end

[var_name,first_,indices] = unique(var_name,'stable');
fprintf(v,'\n$variable names\n');

% do $variable_names next
for ii = 1:length(var_name)
    
    context = var_context{first_(ii)};
    
    if sum(indices == ii) > 1        
      if length(context) > 55, 
           context = sprintf('%s ... (used %d times)',context(1:50), ...
                                                        sum(indices == ii));
      else context = sprintf('%s (used %d times)',context, ...
                                                       sum(indices == ii));
      end
    elseif length(context) > 75, context = [context(1:70) ' ... ']; 
    end
    
    fprintf(v,'\n%s\n%s\n', context, var_name{ii});
end

fprintf(v, '\n\n// End variable definition file\n\n');
% fclose(v);

return



function  [ikey,  ival, dkey] = import_from_variable_file(f_base,varargin)

named = @(v) strncmpi(v,varargin,length(v)); 

do_sweep = any(named('-sweep'));
if do_sweep, ii = find(named('-sweep'));
else         ii = find(named('-usev'));
end

if ii < length(varargin) && exist(varargin{ii+1},'file')
     var_file = varargin{ii+1}; 
else var_file = strrep(f_base,'.template','.variable');
end

v = fopen(var_file,'rt'); 
assert(v > 0, 'failed to open variable file');
active_mode = '#'; 


if any(named('-iter')), iteration = varargin{find(named('-iter'))+1};
  if ischar(iteration), iteration = str2double(iteration); end
else                    iteration = 1; 
end

ikey = {};
dkey = {};
ival = {}; 

if do_sweep, error TODO_run_sweep_from_varfile, end

while ~feof(v)
  
  in = fgetl(v);  
  if ~ischar(in), continue, end
  if isempty(in), continue, end % Blank line   
  in = strtrim(regexprep(in,'//.*$','')); 
  
  if isempty(in), continue, end % Blank line after comment removal
  
  while ~isempty(regexp(strtrim(in),'\.{3}','once')) % "..."
    in = [regexprep(in,'\.{3}.*','') ' ' strtrim(fgetl(v))];
  end
  
  if ismember(in(1),'$#'), active_mode = in(1); continue, end
  if ~any(in == '='), continue, end
  
  key = regexp(in,'^[\w]+','match','once');
  try var = eval(regexp(in,'(?<==)[^%]+','match','once'));
  catch, warning('makeTemplate:varfile:parse','could not evaluate "%s"',in)
    continue
  end
  
  if active_mode == '#' % parsing DEFINE statements
    if var, dkey{end+1} = key; end %#ok<AGROW>
    continue
  else % parsing $variable declarations
    
    if ischar(var), % Use 'string' as is
      ikey{end+1} = key; %#ok<AGROW>
      ival{end+1} = var; %#ok<AGROW>
      continue
    end
    if iteration > numel(var) && numel(var) > 1
        warning('makeTemplate:varfile:itercount', ...
                'possible iteration overflow on iter %d at\n%s',iteration,in)
    end
    var = var(mod(iteration-1,numel(var))+1); % wrap on overflow
    if iscell(var), var = var{1}; end
    ikey{end+1} = key;  %#ok<AGROW>
    ival{end+1} = var;  %#ok<AGROW>
  end
end
fclose(v); 

return
