
function out = cache(mode,newpath) % cache path tools
% This function maintains a local cache in the tempdir or in a specified
% folder for the MATLAB + NEURON pelvic nerve model. A seperate cache is
% maintained for each instance of Matlab; if you close and re-open MATLAB
% you may not be able to find your prior cache
% 
% Usage:
%   cache_folder = tools.cache('PATH') % returns the path to the folder
%   cache_file   = tools.cache('PATH',file_name) % returns specified file    
%                  tools.cache('NEW') % clears the existing cache
%                  tools.cache('NEW',new_path) % make a new empty cache
%                  tools.cache('OLD') % change path to pre-existing cache
%                  tools.cache('SET',new_path) % change path, do not delete
%                                                files in cache
%                  tools.cache('CLEAR',file_spec) % clear selected files
%   file_path    = tools.cache('GET',file_spec) % get specified files
%   file_obj     = tools.cache('GET-STRUCT',file_spec) % as dir() output
%   file_path    = tools.cache('NEWEST',file_spec) % get just newest file
% 
% version 0.3  03-Jun-2020  Calvin Eiber  Refactored for BIDS-like model

if nargin < 1, mode = 'PATH'; end
if nargin < 2, newpath = '*.*'; end

persistent cachepath
if isempty(cachepath)
  
    if isempty(which('getpid')), pid = feature('getpid');
    else                         pid = getpid; 
    end  

    cachepath = 'mdl-VNS-%d'; 
    cachepath = [tempdir sprintf(cachepath,pid)];
    if ~exist(cachepath,'dir'), mkdir(cachepath); end
    
    % if isempty(t),
    % elseif nargin > 1, cachepath = newpath; end
    % else error('%s%s\n\n%s\nparfor ii = 1:N\n\t%s\n\t...\nend', ...
    %          'When using tools.cache inside a parfor loop, ', ...
    %          'you need to explicitly set the cache path e.g:', ...
    %          'cache_path = tools.cache(''reset'');', ...
    %          'tools.cache(''set'',cache_path);')
    % end
end

if ~isempty(which('confirm_recursive_rmdir')) 
  confirm_recursive_rmdir(false); % disable octave confirmation
end

switch(upper(mode))

    case {'PATH','GETPATH'}
        if ~exist(cachepath,'dir'), mkdir(cachepath); end
        if nargin > 1, out = [cachepath newpath]; 
        else           out = cachepath;
        end
        out = regexprep(out,'[/\\]',filesep);
        return        
        
    case {'NEW','RESET'}
        
        if exist(cachepath,'dir'), % remove, using fclose all if needed
          try rmdir(cachepath,'s'); 
          catch C, warning(C.getReport)
            delete(gcp('nocreate')); % shut down parallel pool
            pause(0.1), fclose all; % close all open file streams
            rmdir(cachepath,'s');
          end
        end
        if nargin > 1, cachepath = newpath; end
        if ~ismember(cachepath(end),'/\'), cachepath(end+1) = '/'; end
        if  exist(cachepath,'dir'), rmdir(cachepath,'s'); end 
        if ~exist(cachepath,'dir'), mkdir(cachepath); end
        if nargout > 0, out = cachepath; end
          % out = exist(cachepath,'dir') > 0; end

    case {'SET'} % as NEW but do not delete if exists
        
        if nargin > 1, cachepath = newpath; end
        if ~ismember(cachepath(end),'/\'), cachepath(end+1) = '/'; end
        if ~exist(cachepath,'dir'), mkdir(cachepath); end
        if nargout > 0, out = cachepath; end        
    case 'OLD'
        list = dir(regexprep(cachepath,'[^\\/]+[\\/]$','pn-*'));
        if isempty(list), mkdir(cachepath);
        else sel = (list.datenum == min([list.datenum]));
          cachepath = [list(sel).folder filesep list(sel).name];
        end        
        if nargout > 0, out = cachepath; end
    case {'CLEAR','DELETE'}
        
        if ~exist(cachepath,'dir'), mkdir(cachepath); end
        
        files = dir([cachepath newpath]);
        files([files.isdir]) = []; 
        arrayfun(@(f) delete([cachepath f.name]), files);
        if nargout > 0, out = numel(files); end        
        
    case {'GET','GETFILE'}
        
        files = dir([cachepath newpath]);
        files([files.isdir]) = []; 
        out = arrayfun(@(f) [f.folder filesep f.name], files,'Unif',0);
        if numel(out) == 1, out = out{1}; end

    case {'NEWEST'}
        files = dir([cachepath newpath]);
        files([files.isdir]) = []; 
        files = files([files.datenum] == max([files.datenum])); 
        out = [files(1).folder filesep files(1).name];
        
    case {'GET-S','GET-STRUCT'}
        files = dir([cachepath newpath]);
        files([files.isdir]) = []; 
        out = files;
        
    otherwise warning('PN_Model:cache:unknowncmd', ... 
                      'Unknown request "%s".', mode)
end




