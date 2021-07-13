
function [output] = configuration(varargin)
% open "+tools/configuration.json" and look for an entry matching this
% computer name, username, and machine type. if nargin = 0, returns a 
% structure containing the configuration entries for this machine. If
% nargin > 0, returns the requested entry as a string. 
% (e.g. tools.configuration('eidors') returns the path to the local
%       installation of EIDORS toolbox. string "null" is returned as '')
  
  named = @(v) strncmpi(v,varargin,length(v));
  if nargin > 0, user_key = strrep(varargin{1},'-','_'); end
  
  % if file does not exist, generate
  if any(named('noload')), CONFIG_file = ''; 
  else
    CONFIG_file = tools.file('~/code/+tools/configuration.json');
    if ~exist(CONFIG_file,'file')    
      if isdeployed
          fprintf('+tools.configuration.m: "%s" did not exist, %s ... \n', ...
                         CONFIG_file, 'attempting to generate empty file')
      end
      system(sprintf('touch "%s"',CONFIG_file));
      % fclose(fopen(CONFIG_file,'at'));
    end
    if any(named('open')), edit(CONFIG_file), return, end
  end

  persistent me this
  if any(named('refresh')), me = []; end
  if isempty(me)
    [me,this] = gather_configuration(CONFIG_file); 
  end
  if any(named('make-new-file'))
    CONFIG_file = tools.file('~/code/+tools/configuration.json');
    make_new_config_file(CONFIG_file,me,varargin{:})
    me = []; this = []; return
    
  end
  if any(named('noload')), output = me; return, end
  if isempty(this) 
    [me,this] = gather_configuration(CONFIG_file); 
  end
  
  
  
  if ~isempty(this)
% if found_match % JSON file had a "machine, user, pc-name" match
    if nargin == 0, output = this;
    elseif isfield(this,user_key)
      output = this.(user_key);
      if any(output == '~'), output = tools.file(output); end
      if strcmpi(output,'null'), output = ''; end
    else output = ''; 
    end
    return
  end
  
  warning('ViNERS:config:addingEntry', ...
          'Adding configuration file entry for name="%s", user="%s", machine="%s"', ...
           me.name, me.user, me.machine)
    
  fid = fopen(CONFIG_file,'r+t');
  
  contains = @(a,b) ~isempty(strfind(a,b));
  
  while ~feof(fid)
    p = ftell(fid);
    txt = strtrim(fgetl(fid));
    if ~contains(txt,'}]'), continue, end    
    fseek(fid,p,'bof');
    break
  end
  
  fprintf(fid,'  },{\n');
  for key = {'name','user','machine'}
    fprintf(fid,'    "%s": "%s",\n',key{1},me.(key{1}));
  end
    fprintf(fid,'    "%s": "pn-mdl-%%d/"','cache');
  for key = {'gmsh','neuron','eidors','toolbox-gsa'}
    fprintf(fid,',\n    "%s": "%s"',key{1},'null');
  end  
  fprintf(fid,'  }]\n}');
  fclose(fid);

  if nargin == 0, output = me;
  elseif isfield(me,user_key)
    output = me.(user_key);
    if any(output == '~'), output = tools.file(output); end
    if strcmpi(output,'null'), output = ''; end
  else output = ''; 
  end
  
  clear me this
  return
    
function [me,this] = gather_configuration(CONFIG_file)
  
  if isunix
    me.name = getenv('HOSTNAME');
    me.user = getenv('LOGNAME');    
    me.machine = computer('arch');
  else  
    me.name = getenv('computername');
    me.user = getenv('username');    
    me.machine = computer('arch');
  end
  
  assert(~isempty(me.name) && ~isempty(me.user) && ~isempty(me.machine));
  
  if isempty(CONFIG_file), this = []; return, end
   
  % This is not actually a JSON parser, we just scroll through looking for a 
  % JSON object which matches our system name, username, and machine type
  
  if isempty(strfind(ctfroot, 'MATLAB')) %#ok<*STREMP>
    % contains = @(a,b) cellfun(@(s) ~isempty(strfind(s,b)), a);    
    warning off Octave:regexp-lookbehind-limit
    more off % disable octave paging
  end
  
  this = struct;  
  fid = fopen(CONFIG_file,'rt'); 
  while ~feof(fid)
    
    txt = fgetl(fid);
    if isnumeric(txt), continue, end
    txt = strtrim(txt);
    txt = strtrim(regexprep(txt,'//.*',''));
    if isempty(txt), continue, end
    
    if any(txt == '}') % end-of-object
      
      if ~isfield(this,'name') || ~isfield(this,'user') || ... 
                                  ~isfield(this,'machine'), continue
      end
      
      if strcmpi(this.name,me.name) && strcmpi(this.machine, me.machine) && ...
         strcmpi(this.user,me.user), fclose(fid); return
      end
    end
    
    if  any(txt == '{') % start-of-object
      this = struct; 
    end
  
    if  any(txt == ':') % start-of-field
      name = regexp(txt,'(?<=")[^\"]+(?=":)','match','once');
      value = regexp(txt,'(?<=:\s*")[^\"]*','match','once');      
      if isempty(value), continue, end
      name = strrep(name,'-','_');
      this.(name) = value;
    end
  end
  
  this = []; % return empty handed  
  fclose(fid);
  
  return
  
  
function make_new_config_file(file,me,varargin)

if exist(file,'file')
  warning('ViNERS:config:eraseFile', ...
          'Resetting configuration file. %s for name="%s", user="%s", machine="%s"', ...
          'New file will contain one entry', me.name, me.user, me.machine)
  if isunix
       system(sprintf('mv -f "%s" "%s.save"', file, file));
  else system(sprintf('move /Y "%s" "%s.save"', file, file));
  end
end


is_source = cellfun(@isstruct,varargin); 
if any(is_source), source = varargin{find(is_source,1)}; 
else               source = struct; 
end

fid = fopen(file,'wt'); 

fprintf(fid,'\n// CONFIGURATION file generated %s\n', datestr(now));
fprintf(fid,'\n{ "env-list": [{\n');

for key = {'name','user','machine'}
 fprintf(fid,'    "%s": "%s",\n',key{1},me.(key{1}));
end

[ok,val] = system('git rev-parse --short HEAD'); 
if ok, fprintf(fid,'    "%s": "0.0",\n','version'); % <<<<< UPDATE fixed version here
else   fprintf(fid,'    "%s": "0.0 build %s",\n','version',strtrim(val));
end

fprintf(fid,'    "%s": "pn-mdl-%%d/"','cache');
for key = {'gmsh','neuron','eidors'}    
 if ~isfield(source,key{1}), source.(key{1}) = 'null';end
 fprintf(fid,',\n    "%s": "%s"',key{1},strrep(source.(key{1}),'\','/'));
end  
fprintf(fid,'\n  }]\n}');
fclose(fid);


return




