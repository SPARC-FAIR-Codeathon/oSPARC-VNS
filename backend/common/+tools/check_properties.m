
function check_properties(list,callback)

if ~exist('list','var'), list = tools.file('~/output/'); end
if ~exist('callback','var'), callback = @numel_chunk_labels; end
if ischar(list), list = dir(list); end

isdd_ = @(x) cellfun(@(u) all(u=='.'), {x.name}); 
f_ = @(x,varargin) [x.folder filesep x.name varargin{:}];

list(isdd_(list)) = []; 

files = list(~[list.isdir]); 
list(~[list.isdir]) = []; 

if ~isempty(files), callback(files), end % Check files for what you want to be true 
  
for ii = 1:numel(list)  
  tools.check_properties(dir(f_(list(ii),'/*')))
end


%%


function numel_chunk_labels(files)

f_ = @(x,varargin) [x.folder filesep x.name varargin{:}];

chunk_lbls = unique(regexp({files.name},'chunk=\d+','match','once')); 

if numel(chunk_lbls) == 1    

  files(1).name = chunk_lbls{1}; 
  disp(tools.file('T',f_(files(1))));

elseif numel(chunk_lbls) > 1    
  disp(tools.file('T',f_(files(1))))
  cellfun(@disp,chunk_lbls)
  warning('\nMultiple chunks here')
end