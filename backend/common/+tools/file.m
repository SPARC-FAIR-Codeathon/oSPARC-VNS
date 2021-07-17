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
      warning('BIDSfile:oSPARC_not_supported','SET not supported on oSPARC deployment')
      return
      
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
    case {'ROOT','SET-ROOT'}, update_json(info,'root',varargin{1});
    otherwise
      error('Unknown command "%s"', cmd)
  end
end

%%
function filepath = get_complete_filepath(filename,info,varargin)

named = @(v) strncmpi(v,varargin,length(v));
tok  = find(filename == '~',1);

if filename(1) == '~' 
  root = info(1).root; % e.g. ~/source/***
elseif strncmpi(filename,'out~',4)
  root = [info(1).root filesep 'output'];
elseif strncmpi(filename,'in~',3)
  root = [info(1).root filesep 'input']; 
else
  stub = filename(1:tok-1);  
  for search = {'input','output',''}
    root = [info(1).root filesep search{1} filesep stub]; 
  
    if exist(root,'dir') == 7, break, end
  end
  
  if exist(root,'dir') ~= 7
      error('Neither ''in~\\%s'' nor ''out~\\%s'' exists, %s', ...
              stub, stub, 'please check your requested path')
    
      % going to have to be strict about this, I think.
      error('%s, please refactor "%s" to use in~, out~, or ~\\path', ...
              'This functionality has been removed from tools.file', ...
               filename) %#ok<UNRCH>
  end
end

filepath = [root filename(tok+1:end)]; 
filepath = regexprep(filepath,'[\\/]+',filesep); 

if any(named('prompt')) || any(named('-prompt'))
  filepath = get_filepath_prompt(filepath);
  return
end

if exist(filepath,'file') || exist(filepath,'dir'), return, end
if any(filepath == '*'), return, end % probably going to dir(), do not warn or look for a .ref
if any(filepath == '%'), return, end % probably going to tools.file('GET', ... )

a_dir = ismember(filepath(end),'/\'); % asked for a folder path? 
if a_dir && any(named('-make')), mkdir(filepath), return, end
if a_dir, a_dir = 'folder'; else a_dir = 'file'; end

if any(named('-v')) || isdeployed
    warning('BIDSfile:notFound','Requested %s "%s" does not exist.',a_dir,filename), 
end

return


%% Child functions
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
  if isdeployed, this.root = pwd; 
  else this.root = fileparts(mfilename('fullpath'));  
  end
  
  % this.file = [this.root filesep 'running.json'];
  % this.root = regexprep(this.root,'([\\/])code[\\/].*$','');    
  this.root = regexprep(this.root,'([\\/])common[\\/].*$','');
  if isdeployed
    fprintf('initialising tools.file: ''~'' = ''%s''\n', this.root)
    % fprintf('pwd = %s\nmfilename(''fullphttps://github.com/SPARC-FAIR-Codeathon/oSPARC-VNSath'') = %s\n', pwd, mfilename('fullpath'))
    % fprintf('tools.config(''root'') = %s', tools.configuration('root'))
  end
end

if isempty(strfind(ctfroot, 'MATLAB')) %#ok<STREMP>
  if isdeployed, this.pid = feature('getpid');
  else           this.pid = getpid; 
  end
else this.pid = feature('getpid');
end, this.time = now;

if nargin > 1 % Name, value syntax
  for ii = 1:2:length(varargin)
    this.(varargin{ii}) = varargin{ii+1};
  end
end

% #ifdef OSPARC
info = this;

return
% #else NORMAL_DEPLOY
% ... 135 lines of code removed from tools.file

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

