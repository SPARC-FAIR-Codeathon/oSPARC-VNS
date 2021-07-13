function printInfo(varargin)
% Tools.printInfo(...) operates like fprintf(...), except that it erases the
% previous line of text before writing a new line. 

persistent str

if nargin == 0, str = ''; return, end

if isdeployed, fprintf('\n') 
else, fprintf('%s',char(8*ones(size(str))))
end

str = sprintf(varargin{:});
fprintf('%s',str)