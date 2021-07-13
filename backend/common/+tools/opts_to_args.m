
function list = opts_to_args(args,fields,varargin)
% ViNERS was developed using a consistent 'name',[value] syntax for 
%   controlling the behaviour of most of the modules in +models / +plots. 
% Tools.opts_to_args extends this syntax by enabling model modules to 
%   accept arguments in the form of an options structure. Three equally 
%   valid formats are available:
% 
% %% Example 1: name,value syntax
% models.axon_sfap( '-file', eidors_file, '-axons', axons_file )
% 
% %% Example 2: options structure
% options.sfap.file = eidors_file;
% options.sfap.axons = axons_file;
% models.axon_sfap( options )
% 
% %% Example 3: cell argument syntax
% options.sfap = {'-file', eidors_file, '-axons', axons_file }
% models.axon_sfap( options )
% % equivalent to models.axon_sfap( options.sfap{:} )
% 
% %% Example 4: hybrid syntax
% options = struct; 
% options.sfap.file = eidors_file;
% models.axon_sfap( options, '-axons', axons_file )
% % This will also work with cell argument syntax
% 
% If an argument is supplied in multiple places, the options structure will
% be overridden unless --s2a-first is an argument. --no-s2a will cause
% this to exit after doing nothing. Unless --no-sALL is supplied, the
% options structure 'all' is also included (with the lowest precedence). 
% 
% v0.1 CDE 23-June-2021


if nargin == 0, list = {}; return, end
if isempty(args), list = {}; return, end
if isstruct(args), args = {args}; end
list = args;

named = @(v) strncmpi(v,[args varargin],length(v)); 

if any(named('--no-s2a')), return, end % directive to skip

if ~isstruct(args{1}), return, end
if ~iscell(fields), fields = {fields}; end
if nargin > 2, fields = [fields varargin]; end

if ~any(named('--no-sALL')), fields = [fields {'all'}]; end % directive to skip

opts = args{1}; % input options structure
if ~any(isfield(opts,fields))    
    
  if isfield(opts,'name') && isfield(opts,'folder') && ...
     isfield(opts,'isdir') && isfield(opts,'datenum')
     % return, % s was a the output of dir(); pretty normal. 
  end
    
  fields = sprintf(', %s',fields{:});
  warning('ViNERS:options2arglist', ... 
          'None of "%s" are fields of the input options structure, ignoring.', ...
           fields(3:end))
  return
end

%%
skip_if_found = ~any(named('--s2a-first'));
replace_underscore = ~any(named('--s2a-keep'));

for ff = numel(fields):-1:1 % reverse order
    if ~isfield(opts,fields{ff}), continue, end % missing
    if iscell(opts.(fields{ff})) % {:} expansion 
        list = [list opts.(fields{ff})];  %#ok<AGROW>
        continue
    end
    for name = fieldnames(opts.(fields{ff}))' 
      
      if replace_underscore
        argname = ['-' strrep(name{1},'_','-')];
        argname = strrep(argname,'-[-]+','_');
      else argname = name{1}; 
      end

      argval  = opts.(fields{ff}).(name{1});
      idx = strncmpi(list,argname,numel(name{1})+1);
      if any(idx), 
           if skip_if_found, continue, end
           list(find(idx,1)+1) = opts.(fields{ff}).(name{1});
      elseif ~islogical(argval) || numel(argval) > 1
          list = [list {argname argval}]; %#ok<AGROW>
      elseif argval, list = [list {argname}]; %#ok<AGROW>
      end
    end
end

list(1) = [];

return



    
    
        
        