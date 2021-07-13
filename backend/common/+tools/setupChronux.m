
function opts = setupChronux( varargin )
% Returns some reasonable default options for Chronux toolbox. If Chronux
%   is not on the path, looks for the toolbox in the following: 
%    1. tools.configuration('chronux')
%    2. %userpath%/Add-Ons/Chronux
% 
%   and attempts to add it and sub-folders to the path. 
% 
% Version 0.1  23-Oct-2019  Calvin Eiber & Sander Pietersen

% Default options for Chronux
opts.tapers = [3 5];
opts.Fs = 1000 ;
opts.fpass = [0.5 300];
opts.pad = 0;
opts.err = 0;
opts.trialave = 0;

if evalin('caller','exist(''fs'',''var'')')
    opts.Fs = evalin('caller','fs');
end

if ~isempty(which('mtspectrumc')), return, end

addon_ = @(varargin) fullfile(userpath,'Add-Ons',varargin{:});

% get platform-independent path to LFP analysis toolbox 
p_root = tools.configuration('chronux'); 
if ~exist(p_root,'dir'), p_root = addon_('Chronux'); end
if ~exist(p_root,'dir'), p_root = addon_('CHRONUX_2018'); end
if ~exist(p_root,'dir'), warning('Could not find Chronux'), 
else addpath(genpath(p_root)) % add all, recursively
end






