
function setupEIDORS( varargin )
% If EIDORS is not on the path, looks for the toolbox in the following: 
%    1. %userpath%/Add-Ons/EIDORS
%  and run eidors/startup.m 
% 
% Version 0.2  04-Jul-2020  Calvin Eiber
% - Updated to use tools.configuration(...)

if ~isempty(which('show_fem')),return,end

eidors_path = tools.configuration('eidors'); 

if isempty(eidors_path)
  if ~isempty(which('expdir')), eidors_path = expdir('eidors'); 
  else warning('Missing EXPDIR'), return, 
  end
end
    
if ~exist([eidors_path '/startup.m'],'file')
    warning('eidors/startup.m not found')
    tools.configuration('open');
end

if isempty(which('show_fem')), run([eidors_path '/startup.m']), end
warning('off','EIDORS:normalize_flag_not_set'); % Disable annoying warning