
function make_SPARC_structure(location,varargin)
% Generate core folders (~/source, ~/primary, ~/code, ~/docs, etc)

named = @(v) strncmpi(varargin,v,length(v));
sys_ = @(varargin) system(sprintf(varargin{:}));

if nargin == 0
  location = fileparts(mfilename('fullpath'));
  location = strrep(location,filesep,'/'); % unix-like  
  location = regexprep(location,'/code.*',''); % if /code/..., go to root  
  location = regexprep(location,'/\+.*',''); 
end

here = fileparts(mfilename('fullpath'));
here = strrep(here,filesep,'/'); % unix-like
here = regexprep(here,'/\+tools.*',''); % if /code/..., go to root

dd_ = @(v) cellfun(@(s) all(s=='.'), {v.name}); % ".", ".."
f_ = @(x) [x.folder filesep x.name];

here_list = dir(here); 
here_list(dd_(here_list)) = [];

if ~isfield(here_list, 'folder') % OCTAVE compatiblity
  
  [here_list.folder] = deal(here);
end

% CORE folders for BIDS filestructure
for folders = {'code','source','primary','derivative','docs'}
  p = [location filesep folders{1}]; 
  if ~isfolder(p), mkdir(p); end
  
  here_list(ismember({here_list.name},folders{1})) = []; 
end

here_list = arrayfun(f_,here_list,'unif',0);

if ~any(named('-nomove'))
    % everything 'here' gets moved to /code
    mov_dst = strrep(sprintf('"%s/code/"',location),'/',filesep);
    
    if isunix
         mov_cmd = sprintf('"%s" ', here_list{:});
         sys_('mv %s -t %s',mov_cmd,mov_dst);
    else
      for ii = 1:numel(here_list)
         sys_('MOVE "%s" %s',here_list{ii},mov_dst);
      end
    end
end


%% README.md file

p = [location filesep 'README.md']; 
if ~isfile(p), md = fopen(p,'wt'); 
  w_ = @(varargin) fprintf(md,varargin{:});   
  
  w_('\n# ViNERS dataset \n\n');  
  w_('This document describes the organisaton of `%s` into a SPARC-compatible dataset.\n', location);
  w_('Beware, this does not yet contain SPARC-compliant metadata tables ');
  w_('(subjects.xlsx, samples.xlsx, dataset_description.xlsx, *\\manifest.xlsx)\n');
  w_('Generated %s\nUpdated 2020-Aug-05 Calvin Eiber',datestr(now));
  
  w_('\n\n---\n# DATA folders\n\n## /source/\n');
  w_('Common input to ViNERS\n\n');
  w_('## /primary/\nContains individual model runs, with format\n');
  w_('- `/primary/sub-1/run-xx/axons`: axon populations and membrane currents\n');
  w_('- `/primary/sub-1/run-xx/eidors`: electro-anatomical model results\n');
  w_('- `/primary/sub-1/run-xx/mesh`: input geometry for electro-anatomical model ');
  w_('(optional, defaults to /source/mesh if not found)\n');
  w_('- `/primary/sub-1/run-xx/thresholds`: electrical stimulation trhesolds\n');
  w_('- `/primary/sub-1/run-xx/waves`: simulated nerve recordings \n');
  w_('\n one ''subject'' corresponds to one configuration of fascicles to be simulated. ');
  w_('Within a subject, that arrangement of fascicles can be simulated in multiple ways, ');
  w_('leading to multiple runs. Each electroanatomical ');
  w_('model can be used to generate many recording simulations (e.g. different spike-');
  w_('rates), and so a subject will typically have many files within each sub-/waves/ ');
  w_('folder, organised by simulated population response type.\n\n')  ;  
  w_('## /derivitive/\nContains generated .PDFs from +plots, organised by subject\n\n');
  
  w_('\n\n# SUPPORT folders\n\n## /code/\n');
  w_('Contains the model code, organised into +models, +plots, +utilities, and some example experiments\n\n');
  w_('## /docs/\nContains additional documentation including instructions on how to use the model\n\n');
  w_('## /protocol/ Contains additional information about biological sample collection ');
  w_('and aspects of the experimental design\n\n---\n');
  
  fclose(md); clear w_ md
end



%% Folders in /source for Pelvic Nerve model

% for folders = {'fascicles','axons','mesh'}
%   p = [location filesep 'source' filesep folders{1}]; 
%   if ~isfolder(p), mkdir(p); end
% end

p = [location filesep 'source' filesep 'README.md']; 
if ~isfile(p), md = fopen(p,'wt'); 
  w_ = @(varargin) fprintf(md,varargin{:}); 
  w_('\n# \\source\\\n');
  w_('Optional folder that contains raw data prior to conversion to the ');
  w_('format contained in the primary data folder. We are using this folder ');
  w_('for common input to ViNERS.\n\n');
  
  w_('## \\source\\nerves\\\n');
  w_('splinedata from microscopy sections or generated from the literature ');
  w_('suitable for import into GMSH using mesh.insert_gmsh_fascicles\n');
  
  w_('## \\source\\axons\\\n');
  w_('axon population data from EM or from the literature, needed for +models.axon_population\n\n');
  
  w_('## \\source\\mesh\\\n');
  w_('default generic GMSH .geo model framework, needed for +models.nerve_anatomy.');
  w_('These files get translated to a gmsh .geo file using tools.makeFromTemplate(). '); 
  w_('geometric discriptions of the model features are implemented in +mesh\n');
  w_('Files:\n- pelvic_nerve.geo.template\n- pelvic_nerve.geo.variable\n\n');
  fclose(md); clear w_ md
end

%% Create 'running.json' (multiuser process file)
p = fullfile(location,'code','+tools','running.json'); 
if ~isfile(p) % DNE and needs to be created    
  if isempty(which('tools.configuration'))
    cd([location filesep 'code'])
  end
  this = tools.configuration('noload'); 
    
  md = fopen(p,'wt'); 
  w_ = @(varargin) fprintf(md,varargin{:}); 
  w_('\n// PROCESS file generated %s\n', datestr(now));
  w_('\n{ "name": "%s",\n   "pids":[\n   ', this.name);
  w_('[00000,"sub-1",1000.00101],\n   [00000,"sub-1",1000.00102]\n  ],\n}\n\n');
  fclose(md); clear w_ md
end

%% Add a subject if no subjects exist yet
if ~isfolder([location filesep 'primary' filesep 'Sub-1'])
 if ~isempty(which('tools.make_SPARC_subject')), tools.make_SPARC_subject(1);   
 else mkdir([location filesep 'primary' filesep 'Sub-1'])
 end
end

%% Should do this, but I'm not (#TODO)
warning('Generating subjects/samples/submission.xlsx is left as an exercise for the data wrangler')

return


% for OCTAVE compatability:
function s = isfolder(f), s = exist(f,'dir') == 7; 
function s = isfile(f), s = exist(f,'file') == 2; 