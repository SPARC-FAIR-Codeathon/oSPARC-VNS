function [data, json] = parse_json(json)
% [DATA JSON] = PARSE_JSON(json)
% This function parses a JSON string and returns a cell array with the
% parsed data. JSON objects are converted to structures and JSON arrays are
% converted to cell arrays.
%
% Example:
% google_search = 'http://ajax.googleapis.com/ajax/services/search/web?v=1.0&q=matlab';
% matlab_results = parse_json(urlread(google_search));
% disp(matlab_results{1}.responseData.results{1}.titleNoFormatting)
% disp(matlab_results{1}.responseData.results{1}.visibleUrl)
%
% Downloaded from:
% https://au.mathworks.com/matlabcentral/fileexchange/20565-json-parser
% 
% version 1.2.0.1 Joel Feenstra (MathWorks)
% Copyright (c) 2016, The MathWorks, Inc.
% All rights reserved.

    if exist(json,'file')
        
        if isdeployed, fprintf('Loading %s\n', json), end
        
        fi = fopen(json,'rt');
        fx = onCleanup(@() fclose(fi));
        json = strtrim(fread(fi,inf,'*char')'); 
        json = regexprep(json, ['//[^' newline ']*($|' newline ')'],'');
        clear fx % close file
    elseif any(ismember('{}',json)) % parse json string not a file 
        if isdeployed, fprintf('Parsing "%s ..."\n', json(1:min(end,30))), end
    else error('File not found: %s', json)
    end
    
    data = cell(0,1);
    
    len_j = length(json); 
    
    while ~isempty(json)
        
        [value, json] = parse_value(json);
        data{end+1} = value; %#ok<AGROW>
        
        if len_j <= length(json), 
            error('JSON parse failed near:\n%s',json(1:min(end,50))); 
        end
        len_j = length(json);
    end
    
    data(cellfun(@isempty,data)) = []; 
    if numel(data) == 1, data = data{1}; end
    data = cleanup_json(data);
end
function [value, json] = parse_value(json)
    value = [];
    if ~isempty(json)
        id = strtrim(json(1));
        json(1) = [];        
        json = strtrim(json);
        if isempty(id), return, end
        
        switch lower(id)
            case '"'
                [value, json] = parse_string(json);
                
            case '{'
                [value, json] = parse_object(json);
                
            case '['
                [value, json] = parse_array(json);
                
            case 't'
                value = true;
                if (length(json) >= 3)
                    json(1:3) = [];
                else
                    ME = MException('json:parse_value',['Invalid TRUE identifier: ' id json]);
                    ME.throw;
                end
                
            case 'f'
                value = false;
                if (length(json) >= 4)
                    json(1:4) = [];
                else
                    ME = MException('json:parse_value',['Invalid FALSE identifier: ' id json]);
                    ME.throw;
                end
                
            case 'n'
                value = [];
                if (length(json) >= 3)
                    json(1:3) = [];
                else
                    ME = MException('json:parse_value',['Invalid NULL identifier: ' id json]);
                    ME.throw;
                end
                
            otherwise
                [value, json] = parse_number([id json]); % Need to put the id back on the string
        end
    end
end
function [data, json] = parse_array(json)
    data = cell(0,1);
    while ~isempty(json)
        if strcmp(json(1),']') % Check if the array is closed
            json(1) = [];
            return
        end
        
        [value, json] = parse_value(json);
        
        if isempty(value)
            ME = MException('json:parse_array',['Parsed an empty value: ' json]);
            ME.throw;
        end
        data{end+1} = value; %#ok<AGROW>
        
        while ~isempty(json) && ~isempty(regexp(json(1),'[\s,]','once'))
            json(1) = [];
        end
    end
end
function [data, json] = parse_object(json)
    data = [];
    while ~isempty(json)
        id = json(1);
        json(1) = [];
        
        switch id
            case '"' % Start a name/value pair
                [name, value, remaining_json] = parse_name_value(json);
                if isempty(name)
                    ME = MException('json:parse_object',['Can not have an empty name: ' json]);
                    ME.throw;
                end
                if ~isvarname(name), name = strrep(name,'-','_'); end
                if ~isvarname(name), name = ['p' name]; end %#ok<AGROW>
                data.(name) = value;
                json = remaining_json;
                
            case '}' % End of object, so exit the function
                return
                
            otherwise % Ignore other characters
        end
    end
end
function [name, value, json] = parse_name_value(json)
    name = [];
    value = [];
    if ~isempty(json)
        [name, json] = parse_string(json);
        
        % Skip spaces and the : separator
        while ~isempty(json) && ~isempty(regexp(json(1),'[\s:]','once'))
            json(1) = [];
        end
        [value, json] = parse_value(json);
    end
end
function [string, json] = parse_string(json)
    string = [];
    while ~isempty(json)
        letter = json(1);
        json(1) = [];
        
        switch lower(letter)
            case '\' % Deal with escaped characters
                if ~isempty(json)
                    code = json(1);
                    json(1) = [];
                    switch lower(code)
                        case '"'
                            new_char = '"';
                        case '\'
                            new_char = '\';
                        case '/'
                            new_char = '/';
                        case {'b' 'f' 'n' 'r' 't'}
                            new_char = sprintf('\%c',code);
                        case 'u'
                            if length(json) >= 4
                                new_char = sprintf('\\u%s',json(1:4));
                                json(1:4) = [];
                            end
                        otherwise
                            new_char = [];
                    end
                end
                
            case '"' % Done with the string
                return
                
            otherwise
                new_char = letter;
        end
        % Append the new character
        string = [string new_char]; %#ok<AGROW>
    end
end
function [num, json] = parse_number(json)
  num = [];
  if ~isempty(json)
    % Validate the floating point number using a regular expression
    [s, e] = regexp(json,'^[\w]?[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?[\w]?','once');
    if ~isempty(s)
      num_str = json(s:e);
      json(s:e) = [];
      num = str2double(strtrim(num_str));
    end
  end
end


%% Added CE 13-06-2021


function data = cleanup_json(data)
% convert [[val,val],[val,val]] to matlab matricies (recursively)
    if isstruct(data)
        for f = fieldnames(data)'
            data.(f{1}) = cleanup_json(data.(f{1}));
        end
    elseif iscell(data)        
        for ii = 1:numel(data)            
            data{ii} = cleanup_json(data{ii});
        end
        
        if ~all(cellfun(@isnumeric,data)), return, end        
        s = cellfun(@size,data,'unif',0);
        if numel(unique(cellfun(@numel,s))) > 1, return, end
        s  = cat(1,s{:}); 
        if all(s(:) == 1), data = [data{:}]; return, end
        
        s(:,end+1) = 1; % matrix has dim=1 in the n+1'th dimension
        [~,a] = unique(s,'rows'); 
        if numel(a) > 1, return, end
        c = find(s(1,:) == 1,1);        
        data = cat(c,data{:}); 
    end
end