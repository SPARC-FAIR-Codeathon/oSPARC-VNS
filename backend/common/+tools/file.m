function filepath = file(filename,varargin)
% tools.file is responsible for file-path management and maintaining the 
%   inputs and outputs of ViNERS in a sparc data  structure. 
% 
% tools.file, by itself, returns the path to the installation directory of 
%       ViNERS (e.g. C:\Users\Calvin\Documents\MATLAB\ViNERS). 
% 
% The basic usage is illustrated by >> tools.file('~/source/example.json'),
%   which returns the absolute path to 'example.json' in the source folder.
%   However, 'tools.file' is most useful when working within a subjects
%   (or sample or run) folder, which is illustrated with the following:
% 
% >> tools.file('set','sub-401\sam-1') % set the active subject
% >> tools.file('sub~\axons\axons (1).mat') % get this subject axons file
% 
% The second call to tools.file can be shortened to:  
% >> tools.file('axons~\axons (1).mat') % get axons file for this subject
% 
% tools.file('out~\path\to\object') looks for the subject/sample/run folder
%   as a subfolder of ~/derivative instead of ~/primary. 
% 
% tools.file supports a number of other operations: 
% - tools.file('list', ...)` list subjects
% - tools.file('list', '<~ path>') list contents of folder
% - tools.file('rename','<~ path>') rename file or folder
% - tools.file('info') return info about current file/subject and
%                        machine/process information
% - tools.file('get','<~ path>') returns a full filepath from a partial 
%                                  filepath. In the example above, 
%                                  tools.file('get','axons~\axons*.mat')
%                                  would return the newest axons.mat file 
%                                  in axons~/ from the current subject. 
% - tools.file('get','<~ path>','next') returns a full filepath for 
%                                  creating a new file which will not 
%                                  overwrite any other files. 
% - tools.file('cache') calls tools.cache, and is provided for convenience
% - tools.file('open','<~ path>') calls winopen on '<~ path>'
% - tools.file('short','path') returns the <~ path> for a given full path, 
%                    with ViNERS' installation directory replaced with '~'
% - tools.file('shorter','path') returns a shorter <~ path> for a given 
%                  full path, with tools.file('sub~') replaced with 'sub~'
% 
% Last updated CDE 22-June-2021


persistent info, info = update_json(info);
if nargin == 0, filepath = info(1).root; return, end
if contains(filename,filesep) && exist(filename,'file') == 2, return, end 

%% Pattern 1 : ~/source/path/to/goods.csv
if any(filename == '~')  
  filepath = get_complete_filepath(filename,info,varargin{:});   
else
  %% Other Instruction
  cmd = upper(filename); 
  switch cmd

    case {'S','SET'}

      if iscell(varargin{1}), varargin{1} = varargin{1}{1}; end
      if isstruct(varargin{1}) && isfield(varargin{1},'folder')
          sub = regexp(varargin{1}.folder,'(?<=[\\/])sub-[^\\/]*','match','once');
          sam = regexp(varargin{1}.folder,'(?<=[\\/])sam-[^\\/]*','match','once');
          run = regexp(varargin{1}.folder,'(?<=[\\/])run-[^\\/]*','match','once');
          varargin = [{sub},{sam},{run},varargin(2:end)];
      end
        
      sub = set_subject_ID(info,varargin{:});
      info = update_json(info,'sub',sub);
      if nargout > 0, filepath = tools.file('sub~/','-make'); end % make if dne

    case {'L','LIST'}, list_path_contents(varargin{:}); 
    case {'H','HELP'}, doc(mfilename) 
    case {'R','RN','REN','RENAME'}, rename_file(varargin{:}); 
    case {'I','INFO'}, filepath = info;
    case {'G','GET'}, filename = varargin{1}; 
      if any(filename == '~'), filename = tools.file(filename,varargin{2:end}); end
      filepath = get_partial_file(filename,varargin{2:end}); 
    case {'C','CACHE'}, filepath = tools.cache(varargin{2:end});
    case {'O','OPEN','GOTO'}, winopen(tools.file(varargin{:})); 
    case {'T','TITLE','SHORT'}, filepath = strrep(varargin{1},info.root,'~');        
    case {'T2','TT','SHORTER'}, filepath = strrep(varargin{1},tools.file('sub~'),'sub~');
    case {'MAKE'}, tools.make_SPARC_subject(varargin{:}); 
                   tools.file('set',varargin{:}); 
    otherwise
      error('Unknown command "%s"', cmd)
  end
end

%%
function filepath = get_complete_filepath(filename,info,varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
if any(named('sub')) || any(named('tag'))
  info = set_subject_ID(info,varargin{:});
  info = update_json(info,'sub',sub);  
end

tok  = find(filename == '~',1);

if filename(1) == '~' 
  root = [info(1).root]; % ~/source/***    
elseif strncmpi(filename,'sub~',4) % ~/primary/sub-XX/folder/data.mat  
  root = [info(1).root filesep 'primary' filesep info.sub];
elseif strncmpi(filename,'out~',4) || strncmpi(filename,'deri',4)
  root = [info(1).root filesep 'derivative']; 
else % other folders inside of 
  
  root = filename(1:tok-1);
  list = dir([info.root '/primary/' info.sub '/*']); 
  list = list([list.isdir]);  
  
  if ~any(strcmp({list.name},root)) && ~any(named('-make'))
    % if not found in /primary/sub-xxx, look in ~/
    list = dir([info.root '/*']); 
    list = list([list.isdir]);
  end
  
  if ~any(strcmp({list.name},root))
    if any(named('-make')) 
      mkdir([list(1).folder filesep root])
      list(1).name = root; 
    else
     if ~any(named('-q'))
       warning('ViNERS:file:missingSampleType', ...
               '"%s" not found in ~/primary/%s or ~/', root, info.sub);
     end
     filepath = ''; 
     return
    end
  end
  root = [list(1).folder filesep root];
end

filepath = [root filename(tok+1:end)]; 
filepath = regexprep(filepath,'[\\/]',filesep); 

if any(named('prompt')) || any(named('-prompt'))
  filepath = get_filepath_prompt(filepath);
  return
end

if exist(filepath,'file') || exist(filepath,'dir'), return, end
if any(filepath == '*'), return, end % probably going to dir(), do not warn or look for a .ref
if any(filepath == '%'), return, end % probably going to tools.file('GET', ... )

a_dir = ismember(filepath(end),'/\'); % asked for a folder path? 
if a_dir && any(named('-make')), mkdir(filepath), return, end

% 
list = dir([fileparts(filepath) filesep '.reference']); 
if ~isempty(list)
  list = list(max_(list,'datenum'));
  fid = fopen([list.folder filesep list.name],'rt');
  lnk = ''; 
  while ~feof(fid)
    lnk = strtrim(fgetl(fid));
    if isempty(lnk), continue, end
    if ismember(lnk(1),'%#'), lnk = ''; continue, end        
    if any(lnk == '~'), lnk = tools.file(lnk); end
    break
  end
  fclose(fid);
  if isempty(lnk)
      warning('ViNERS:file:emptyRefFile', ...
          'The path "%s" contained a .reference file, %s', p_path, ...
                    'but that .reference file was empty!')
  elseif ~exist(lnk,'dir')
      warning('ViNERS:file:badRefFile', ...
          'The path "%s" %s, but the linked path "%s" was not found.', ...
                    p_path, 'contained a .reference file', lnk);
  else filepath = strrep(filepath,fileparts(filepath),lnk);
  end
end

if a_dir, a_dir = 'folder'; else a_dir = 'file'; end
if any(named('-v')), warning('ViNERS:file:notFound','Requested %s "%s" does not exist.',a_dir,filename), end
return


%% Child functions
function [name,list] = set_subject_ID(info,varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
name = info(1).sub; 

% if ~any(varargin{1} == '='), varargin{1} = strrep(varargin{1},'-','='); end

list = dir([info(1).root '/primary/sub*']); % try lowercase first
if numel(list) == 0, list = dir([info(1).root '/primary/Sub*']); end
if numel(list) == 0, list = dir([info(1).root '/primary/pool*']); end
if numel(list) == 0, list = dir([info(1).root '/primary/*']); end
  
list(~[list.isdir]) = []; % directories only
max_ = @(f,v) [f.(v)] == max([f.(v)]);

if isempty(list)
  error('could not find anything ~/primary. Try running %s,%s.', ...
        'tools.make_SPARC_structure','tools.make_SPARC_subject(1)')
end

[list.id] = deal(0); 
list(1).id = cellfun(@(n) str2double(regexp(n,'\d+','match','once')), ...
                                           {list.name},'unif',0);
[list.id] = list(1).id{:};


if nargin <= 1 % not specified, give a menu and warn that this is not OK
  if all(cellfun(@(x) numel(x)==1, {list.id}))
    [~,seq] = sort([list.id]); list = list(seq); 
  end  
  
  warning('ViNERS:file:pleaseSetSubject', ...
          'Please call tools.file(''set'',subject-id)')
  sel = menu('Select active subject',{list.name});
  if sel > 0, name = list(sel).name; end
  return
end

%%

if ~isnumeric(varargin{1}) % parse string argument 
  
  cmd = varargin{1}; 
  if any(ismember(cmd,'/\')) % filepath fragment    
    cmd = strsplit(cmd,{'\','/'});    
    if any(strncmpi(cmd,'run',3))
        varargin = [cmd(strncmpi(cmd,'run',3)) varargin]; 
    end
    if any(strncmpi(cmd,'sam',3))
        varargin = [cmd(strncmpi(cmd,'sam',3)) varargin]; 
    end 
    varargin = [cmd(1) varargin];    
    named = @(v) strncmpi(v,varargin,length(v)); % update, varargin changed
    cmd = cmd{1};
  end
  cmd = strtrim(regexp(cmd,'(?<=[=\-:]).*','match','once'));
  if ~isnan(str2double(cmd)), varargin{1} = str2double(cmd);
  elseif ~isempty(dir([info.root '/primary/' varargin{1} '*'])) 
    % exact or prefix match    
    list = dir([info.root '/primary/' varargin{1} '*']); 
    [~,sel] = min(cellfun(@length,{list.name}));
    if numel(list) ~= 1
      warning('ViNERS:file:multipleSubjectsListed', ...
              '%d subjects in ~/primary/ match "%s"',numel(list),varargin{1})
    end
    name = strrep([list(sel).folder filesep list(sel).name],  ...
                  [info.root '/primary/'],'');

    if any(named('sam')) % No error/existance checking
      sam = varargin{named('sam')}; 
      if any(sam=='='), sam = regexp(sam,'(?<=[=]).*','match','once'); end
      name = [name filesep sam]; 
    end
    if any(named('run')) % No error/existance checking
      run = varargin{named('run')}; 
      if any(run=='='), run = regexp(run,'(?<=[=]).*','match','once'); end
      name = [name filesep run]; 
    end
    name = strrep(tools.file('T',name),['~' filesep 'primary' filesep],'');
    return
  elseif any(named('-make'))
      
    name = varargin{1};
    if any(named('sam')) % No error/existance checking
      sam = varargin{named('sam')}; 
      if any(sam=='='), sam = regexp(sam,'(?<=[=]).*','match','once'); end
      name = [name filesep sam];       
    end    
    if any(named('run')) % No error/existance checking
      run = varargin{named('run')}; 
      if any(run=='='), run = regexp(run,'(?<=[=]).*','match','once'); end
      name = [name filesep run];       
    end    
    mkdir([info.root '/primary/' name])

  end
  if isnan(str2double(cmd))
   if iscell(cmd), cmd = cmd{1}; end
   switch lower(cmd(1:min(4,end)))
    case 'curr' % pass
    case 'newe', name = list(max_(list,'id')).name;
    case 'new'
        error('Please call tools.make_SPARC_subject')
    otherwise
      warning('ViNERS:file:unableToParse', ...
          'Unable to parse "%s", ignoring ...', varargin{named('sub')})
   end
  end
end

if isnumeric(varargin{1})
  
  sub = varargin{1};   
  sel = ([list.id] == floor(sub)); 
  if ~any(sel)
    if any(named('-make'))  
      mkdir(fullfile(info.root,'primary',sprintf('sub-%d',floor(sub))));
      tools.make_SPARC_subject(floor(sub));
      [name,list] = set_subject_ID(info,varargin{:}); % try again
      return
    else
      error('Missing sub-%d in ~/primary, %s', floor(sub), ...
            'did you mean to call tools.make_SPARC_subject?')
    end
  end
  
  if sum(sel) > 1
    len = cellfun(@length,{list(sel).name});
    sel(sel) = (len == min(len));
    sel = find(sel,1);
        
    warning('ViNERS:file:multipleSubjectIds', ...
            'Multiple directories in ~/primary match #%d, using "%s"', ...
               floor(sub), list(sel).name)
  end
  name = list(sel).name;
  
  if any(named('sam')) % No error/existance checking
    sam = varargin{named('sam')}; 
    if any(sam=='='), sam = regexp(sam,'(?<=[=]).*','match','once'); end
    name = [name filesep sam];     
  end
  
  if any(named('run')) % No error/existance checking
    run = varargin{named('run')}; 
    if any(run=='='), run = regexp(run,'(?<=[=]).*','match','once'); end
    name = [name filesep run];       
  end
  
  if mod(varargin{1},1) == 0, return, end 
  
  %% Repeat to get .sample for subject.sample syntax
  
  sam = mod(sub,1); 
  if nargin > 2 && isnumeric(varargin{2}), sam = sam * 10.^varargin{2}; 
  else while abs(sam-round(sam)) > 1e-12, sam = 10*sam; end
       sam = round(sam);
  end
  
  list = dir([info.root '/primary/' name '/*']);
  list(~[list.isdir]) = []; % directories only

  [list.id] = deal(0);
  list(1).id = cellfun(@(n) str2double(regexp(n,'\d+','match','once')), ...
                                         {list.name},'unif',0);
  [list.id] = list(1).id{:};

  sel = ([list.id] == floor(sam)); 
  if ~any(sel)
    if any(named('-make')) 
      mkdir([list(1).folder filesep sprintf('sam-%d',sam)])
      
      [name,list] = set_subject_ID(info,varargin{:}); % try again
      return
    else
      warning('ViNERS:file:missingSampleNo', ...
              'Missing sam-%d in %s.', sam, name)
    end
  end

  if sum(sel) > 1
    len = cellfun(@length,{list(sel).name});
    sel(sel) = (len == min(len));    
    sel = find(sel,1); 
    warning('ViNERS:file:multipleSampleIds', ...
            'Multiple directories in %s match #%d, using "%s"', ...
               name, sam, list(sel).name)
  end
  
  if any(sel), name = [name filesep list(sel).name]; end
end
    
if isempty(name)
  name = list(max_(list,'id')).name;
  warning('ViNERS:file:subjectNotSet', ...
          'subject ID was not set, using "%s"', name)
end
return


function fpath = get_filepath_prompt(fpath)

named = evalin('caller','named'); 
a_dir = ismember(fpath(end),'/\'); 

if any(named('-prompt-m')), opt = {'Multiselect','on'};
else opt = {}; 
end

if a_dir, fpath = uigetdir(fpath,opt{:});        
   if all(fpath == 0), fpath = ''; return, end  
else
  [s_path,s_name,s_ext] = fileparts(fpath);    
  [s_name,s_path] = uigetfile(['*' s_ext],'', ... 
                              [s_path filesep s_name s_ext],opt{:});
  if all(s_path == 0), fpath = ''; return, end
  fpath = strcat(s_path,s_name);
end

function info = update_json(info,varargin)

persistent this
named = @(v) strncmpi(v,varargin,length(v)); 

if any(named('clear')), this = []; end
if isempty(this)  
  this = tools.configuration('noload'); 
  this.root = fileparts(mfilename('fullpath'));  
  this.file = [this.root filesep 'running.json'];
  this.root = regexprep(this.root,'([\\/])code[\\/].*$','');  
end
if isempty(strfind(ctfroot, 'MATLAB')), this.pid = getpid; %#ok<STREMP>
else  this.pid = feature('getpid');
end,  this.time = now;

FLAG_update = 0; 
json = [];  %#ok<NASGU>

try tools.semaphore('wait',LOCK_HANDLE);
    tools.semaphore('post',LOCK_HANDLE); 
catch C, tools.semaphore('create',LOCK_HANDLE,1); %#ok<NASGU>
end

if isempty(info) % Persistent info not supplied
  
  if exist(this.file,'file')
    json = json_to_struct(this.file);
    % No longer includes a name check
    if any(json.pid == this.pid) % 
        this.sub = json.sub{json.pid == this.pid};
    else [~,sel] = max(json.time); % most recent
      this.sub = json.sub{sel};        
      json.pid(end+1) = this.pid;
      json.sub{end+1} = this.sub;
      json.time(end+1) = this.time;
      FLAG_update = 1;
    end
  else 
    warning('ViNERS:file:missingProcessIDFile', ...
          '[%d] Could not locate (%s)\n', this.pid, ...
                         strrep(this.file,this.root,'~'))
    this.sub = 0;
  end
else % Persistent info supplied, double-check PIDs
  if info.pid == this.pid 
    if this.time > info.time + 1, FLAG_update = 1; end    
  else FLAG_update = 1;
    error TODO_resolve_check_PID
  end
end

assert(~isempty(this),'"this" became empty after loading process info')

%% Do any updating? 
if nargin > 1 % Name, value syntax
  for ii = 1:2:length(varargin)
    this.(varargin{ii}) = varargin{ii+1};
  end  
  FLAG_update = 1;
end

info = this; % value to be returned
if ~FLAG_update && ~isempty(fieldnames(info)), return, end

tools.semaphore('wait',LOCK_HANDLE);
lock = onCleanup(@()tools.semaphore('post',LOCK_HANDLE));

if ~exist(this.file,'file'), make_pids_file(this.file, info), end

fid = fopen(this.file,'r+t');
str = fread(fid,inf,'*char')';  
json = json_to_struct(str,true);
pad_len = ftell(fid);
frewind(fid);

%% Write JSON file 

timed_out = (json.time < now-7);
json.pid(timed_out) = [];
json.sub(timed_out) = []; 
json.time(timed_out) = [];

this_session = (json.pid == this.pid);
if any(this_session)
  json.sub(this_session) = {this.sub};
  json.time(this_session) = this.time;
else
  json.pid(end+1) = this.pid;
  json.sub(end+1) = {this.sub};
  json.time(end+1) = this.time;
end

fprintf(fid,'\n// PROCESS file generated %s by %s\n', ...
                         datestr(this.time),this.user);
fprintf(fid,'{ "name": "%s",\n  "pids":[\n', this.name);
for ii = 1:numel(json.pid)
  
  if ii < numel(json.pid), comma = ','; else comma = ''; end
  fprintf(fid,'\t[%d,"%s",%f]%s\n',json.pid(ii), json.sub{ii}, ...
          json.time(ii),comma);
end
fprintf(fid,'  ],\n}\n');

dat_len = ftell(fid);
if dat_len < pad_len % add blank to end of file
  fprintf(fid,'%s',char(' '*ones(1,pad_len-dat_len)));
end

clear lock % tools.semaphore('post',LOCK_HANDLE);

return

function s = json_to_struct(filename,recurse_flag)
%%
if nargin == 1 || ~recurse_flag
  
  tools.semaphore('wait',LOCK_HANDLE);
  lock = onCleanup(@()tools.semaphore('post',LOCK_HANDLE));
  
  fid = fopen(filename,'rt');
  json = fread(fid,inf,'*char')';  
  fclose(fid);
  clear lock % unlock
  
else json = filename; 
end

json = regexp(json,'{.*}','match','once'); % remove header footer
json = regexprep(json,'[\\/]+','\\\\');
json = tools.parse_json(json);

json = cat(1,json.pids{:});

s.pid = cat(1,json{:,1});
s.sub = cat(1,json(:,2));
s.time = cat(1,json{:,3});

return


function make_pids_file(filename, info)

fi = fopen(filename,'wt');
fprintf(fi,'\n// PROCESS file generated %s by %s\n',datestr(info.time),info.user);
fprintf(fi,'{ "name": "%s",\n  "pids":[\n', info.name);
fprintf(fi,'    [%d,"%s",%f],\n[%d,"%s",%f]\n  ],\n}', ...
                0,info.sub,info.time,info.pid,info.sub,info.time);
fclose(fi);


function k = LOCK_HANDLE, k = 905792; % randomly generated number

function list = list_path_contents(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 

if nargin > 0 && any(varargin{1} == '~')
     list = dir(tools.file(varargin{1}));
else list = dir(tools.file('~/primary/'));
end

list(cellfun(@(x) all(x=='.'),{list.name})) = [];    
f_ = @(x) [x.folder filesep x.name];

if all(strncmpi({list.name},'sub-',4))
  val = str2double(regexp({list.name},'\d+','match','once'));
  [~,order] = sort(val); 
  list = list(order); 
end

if any(~[list.isdir])
  fprintf('./:%s\n',sprintf('  %s',list(~[list.isdir]).name))
end

for ff = 1:length(list)
  if ~list(ff).isdir, continue, end  
  inner = dir(f_(list(ff)));
  inner(cellfun(@(x) all(x=='.'),{inner.name})) = []; 
  if ~any(named('sub~')), inner(~[inner.isdir]) = []; end

  inner = {inner.name}; fix = contains(inner,' ');
  inner(fix) = cellfun(@(x) ['"' x '"'], inner(fix),'Unif',0);      
  fprintf('%s:%s\n',list(ff).name,sprintf('  %s',inner{:}))
end

if nargout == 0, clear, end

function rename_file(folder,new_name)% from "renameLastFile.m"

if nargin == 1, new_name = folder; folder = pwd; end
  
ext = regexp(new_name,'\..*$','match','once');
n = 1;
while exist(new_name,'file')
    n = n+1;
    new_name = regexprep(new_name,'(\(\d+\))?\.',['(' num2str(n) ').']);
end

files = dir([folder filesep '*' ext]);
[~,i] = max(cellfun(@datenum,{files.date}));
system(['ren "' files(i).name '" "' new_name '"']);

%% from tools.getFile
function filename = get_partial_file(filename, mode, varargin)

if nargin == 1, mode = 'newest'; end
if strncmpi(mode,'-q',2), mode = 'newest'; varargin(end+1) = {'-q'}; end
do_verbose = ~any(strncmpi(varargin,'-q',2));

if ~isempty(varargin) 
  filename = regexprep(filename,'(?<=[^\\])\\(?=[^\\])','\\\\');
  filename = sprintf(filename,varargin{:}); 
end

if strcmpi(mode,'next'), filename = get_partial_next(filename);
  return
end


d = dir(filename); 

fp = fileparts(filename); 
if isempty(fp), fp = '.'; end
if isempty(d) 
  if do_verbose
    warning('getfile:fileNotFound','Could not find a file "%s"', filename), 
  end
  filename = []; 
  return
end

switch lower(mode(1:3))
  case 'new', [~,sel] = max([d.datenum]);     
  case 'old', [~,sel] = min([d.datenum]); 
  case 'big', [~,sel] = max([d.bytes]); 
  case 'sma', [~,sel] = min([d.bytes]); 
  otherwise    sel = 1;
    if do_verbose, warning('getfile:unknownMode','unknown mode "%s"', mode)
    end    
end
  
filename = [fp filesep d(sel).name];  
function filename = get_partial_next(filename)

[fp,fn,fe] = fileparts(filename); 
if isempty(fp), fp = '.'; end
fp = [fp filesep];

if contains(fn,'%t')
  fn = strrep(fn,'%t',datestr(now,'yyyymmdd@HHMMSS'));
  if isempty(dir([fp regexprep(fn,'@[\d]+([\)\]])?','$1') fe]))
       fn = regexprep(fn,'@[\d]+([\)\]])?','$1'); % remove (HHMMSS) if not needed
  else fn = regexprep(fn,'@(\d+)','-$1');
  end
elseif contains(fn,'%d')
  
  fn_star = regexprep(fn,'\s*\(?%d\)?','*'); 
  d = dir([fp fn_star fe]); 
  
  if isempty(d) % fn = regexprep(fn,'[\s]+[^\s]*%d[^\s]*\',''); 
      fn = strrep(fn,'%d',num2str(1));
  else
    fn_star = regexprep(fn,'.*(.)%d(.).*','''$1(\\d+)''$2');
    fn_star = regexprep(fn_star,'''([\[\(\)\]\\])','\\$1');    
    val = regexp({d.name},fn_star,'tokens','once'); 
    val = cellfun(@(v) str2double(v{1}),val(~cellfun(@isempty,val)));
    val = nanmax(val) + 1; 
    if isempty(val) || isnan(val), val = 1; end    
    fn = strrep(fn,'%d',num2str(val));
  end
elseif exist(filename,'file')
  warning('nextfile:overwrite','Preparing to overwrite %s',filename)
end

filename = [fp fn fe];

