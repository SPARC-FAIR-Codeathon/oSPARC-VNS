
function output = INPUT_file(list,object,varargin)

if nargin == 0, object = tools.file('sub~'); end

if nargin == 1, 
  object = list; list = dir; list(:) = [];  
  varargin = {'-nosave'};
end

named = @(v) strncmpi(v,varargin,length(v));
f_ = @(x) [x.folder filesep x.name];

%% Find INPUTS file as close to the level specified, writing one in /sub~/ if none found

IN_file = '';
[o_path,o_lvl] = fileparts(object); 
if isempty(o_path), o_path = tools.file('sub~/'); end

while isempty(IN_file)
  
  if exist([o_path filesep 'INPUTS'],'file')
    IN_file = [o_path filesep 'INPUTS'];
    break
  end
  
  if strncmpi(o_lvl,'sub-x',4)
    IN_file = [o_path filesep o_lvl filesep 'INPUTS'];
    fclose(fopen(IN_file,'at')); % make if it doesn't exist
    break
  end
  
  [o_path,o_lvl] = fileparts(o_path); 
end

if nargin == 0, output = IN_file; return, end

%%

if numel(list) == 1, output = f_(list);
  if any(named('-nosave')), return, end
else output = ''; 
end


fid = fopen(IN_file,'rt');
IN_text = fread(fid,inf,'*char')';
fclose(fid);

IN_text = strsplit(IN_text,newline);
IN_trim = strtrim(IN_text); 
 

p_root = tools.file('sub~\'); 
o_short = strrep(object,p_root,'sub~\');

line = strncmpi(IN_trim,o_short,numel(o_short));

if isempty(output), % 
 if any(line), % We found an entry! use it. (overrides singleton?)
  
  if sum(line) > 1
    error TODDO_add_logic_here_to_pick_relevent_item_multiple_inputs
  end
    
  line = strsplit(IN_trim{find(line,1)},'\t');
  if any(line{2} == '~'), line{2} = tools.file(line{2}); end
  if exist(line{2},'file') == 2, output = line{2};
  elseif ~isempty(list) && exist([list(1).folder filesep line{2}],'file')
    output = [list(1).folder filesep line{2}];
  end
  return
 elseif isempty(list), output = ''; return
 else
  sel = menu(sprintf('Choose file to associate with\n%s',o_short), {list.name});
  output = f_(list(sel)); 

  if any(named('-nosave')), return, end
  IN_text{end+1} = sprintf('%s\t%s\t# added by user, %s',o_short,list(sel).name,datestr(now));
 end
elseif ~any(line), % (list) was a singleton, write new line here for future reference
  IN_text{end+1} = sprintf('%s\t%s\t# auto-gen, %s',o_short,list(1).name,datestr(now)); 
end

% Write new information to file so you can find it next time
fid = fopen(IN_file,'wt');
fprintf(fid,'%s\n',IN_text{:});
fclose(fid);


end