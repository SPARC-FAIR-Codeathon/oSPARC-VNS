
function list = gather(source,pattern,varargin)
% list = tools.gather(source,pattern) 
% Recursively search the folders specified by SOURCE for files matching
%   PATTERN; return these files as LIST. if PATTERN is a function handle,
%   recursively apply that function handle to the files in SOURCE and every
%   subdirectory of SOURCE (no direct output will be generated)
% 
% Source may be: a TABLE containing 'subject' (and optionally 'sample')
%                a STRUCT containing 'name' and 'folder' (see dir)
%                a STRING path or tools.file path, input to dir()
% 
% V0.1 CDE 26-Jan-2021


if nargin < 1, source = dir(tools.file('primary~')); end  
if isa(source,'table') source = table2struct(source); end
if ischar(source), 
  if any(source == '~'), source = tools.file(source); end
  source = dir(source); 
end
if nargin < 2, pattern = 'sens*.mat'; end

if isa(pattern,'function_handle'), search(source,pattern), return, end

locate('RESET'); 

if isfield(source,'name') && isfield(source,'folder')
  tools.check_properties(source,@(p) locate(p,pattern) )
  clc
elseif isfield(source,'subject') && isfield(source,'sample')
  for ii = 1:length(source)    
    tools.file('set',source(ii).subject, source(ii).sample)
    search(tools.file('sub~'),@(p) locate(p,pattern) )
  end
elseif isfield(source,'subject')
  for ii = 1:length(source)    
    tools.file('set',source(ii).subject)
    search(tools.file('sub~'),@(p) locate(p,pattern) )
  end
end

list = locate('GET');

u = arrayfun(@(f) [f.folder f.name], list,'unif',0);
[~,u,~] = unique(u);

list = list(u); 


function search(list,callback)

if ischar(list), 
  if any(list == '~'), list = tools.file(list); end
  list = dir(list); 
end

isdd_ = @(x) cellfun(@(u) all(u=='.'), {x.name}); 
f_ = @(x,varargin) [x.folder filesep x.name varargin{:}];

list(isdd_(list)) = []; % remove '.', '..'

files = list(~[list.isdir]); % files get callback applied to them
list(~[list.isdir]) = [];    % remove files from list 

if ~isempty(files), callback(files), end % Check files using callback 
  
for ii = 1:numel(list) % recursive search 
  search(dir(f_(list(ii),'/*')),callback)
end


function result = locate(list,pattern,varargin)

persistent r

if ischar(list) && nargin == 1
  if strcmpi(list,'RESET'), r = []; return, end
  if strcmpi(list,'GET'), result = r; return, end
end

% if isempty(regexp(list(1).folder,...
%            'primary[/\\]sam.*[/\\]sub[^/\\]+$','once'))
%   return
% end

these = dir([list(1).folder filesep pattern]); 
if isempty(these), return, end

r = [r; these];



  
  
  
  
  
  