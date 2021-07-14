

function EIDORS_fwd_model (mesh_file, sigma_file, varargin)

if nargin < 1 || isempty(mesh_file)
  fprintf('Arg 1 not set, using example mesh file\n')
  mesh_file = './input/demo/example cuff (2).mat'; 
end

if nargin < 1 || isempty(mesh_file)
  % fprintf('Arg 2 not set, using default conductivities\n')
  sigma_file = './input/default-conductivity.json'; 
end

tools.file('root',pwd); % set 'root' to this folder

sigma = tools.parse_json(sigma_file); 
sigma = [sigma.sigma{:}];

if isdeployed, disp('Progress: 10%'), end

run_FWD_model(mesh_file, '-sigma', sigma, varargin{:})


return

function run_FWD_model(varargin)
% EIDORS_fwd_model(geometry, ...) runs the electroanatomical model implemented
% using EIDORS for electrical stimulation & recording from a peripheral 
% nerve, including meshing the model in gmesh and gnerating any neccessary
% thin-layers
% 
% This is iteration #3 of this code, which has fixed issues with lost
% mesh points and uses ./+mesh/pelvic_nerve.geo.template
% (future versions will have more ways of customising the models). 
% 
% "geometry" is a struct which defines the stimulating / recording
% electrode geometry, some meshing parameters, and which fascicle set to
% use (from the default splines.dat, see read_dat_file). If "geometry" is
% not supplied, it's read from tools.cache('PATH','elec-geom.mat')
% 
% Some additional -options can be specified: 
% -regenerate: if set, the mesh is regenerated even if
%              tools.cache('path')/pelvic_nerve-thin.msh.mat exists. 
% -mesh:     generate mesh only and exit without running EIDORS. 
% -fascicle: which fascicles from file get simulated (default: 1-5) 
% -stimulus: compute the electrode-tissue relationship for electrical 
%            stimuli (if not set, this code computes recording sensitivity)
%     -mono: if set, stimuli are monopolar (default: sequential bipoles)
% -reference: default [ COUNTER = 1mm, -1, d=2   GROUND = 1mm, -1, d=3 ] 
% -common: use the COUNTER for both COUNTER and GROUND
% -output: sets output file, default '~/eidors/sensitivity (%d).mat'
% 
% Calvin Eiber 25-Apr-2020

%% Parse model inputs

tools.setupEIDORS;

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

do_PLOT = any(named('-p')) == isdeployed; 

if any(named('-q')) && ~isdeployed, printf = @(varargin) []; 
else                                printf = @fprintf; 
end

%% Load specififed mesh

if any(named('-mesh')), input_mesh = get_('-mesh');
elseif nargin > 0 && exist(varargin{1},'file'), input_mesh = varargin{1};
else error('No mesh file supplied, please run module_nerve_mesher')
end
    
m = load(input_mesh); 

%% Construct virtual point electrodes and stimulation pattern for computation of sensitivity

printf('Constructing current pattern ... \n')
% note: this code depends on exact m.object_names 
if isdeployed, disp('Progress: 15%'), end

m.electrode(cellfun(@isempty,{m.electrode.name})) = [];
nE = numel(m.electrode);

if do_PLOT
    clf, show_fem(m), hold on
    h = get(gca,'Children'); 
    set(h,'EdgeAlpha',0.3)
    set(h,'FaceAlpha',0.1)
end

if any(named('-set-zc'))
  [m.electrode.z_contact] = deal(get_('-set-zc'));
end

v_ = @(x) reshape(x,1,[]);

%% Set stimulating and reference electrodes

ref = nE + 1;
meas_pattern = sparse(v_((1:nE)' * [1 1]), ...
                         [1:nE ref*ones(1,nE)], ...
                         [ones(1,nE) -ones(1,nE)]); % Updated reference
m.stimulation = struct([]);

for ee = 1:nE 
    m.stimulation(ee).stimulation = 'Amp';
    stim_pattern = sparse([ee ref],[1 1],[1 -1],nE+1,1); % Monopolar pattern 
    m.stimulation(ee).stim_pattern = stim_pattern;
    m.stimulation(ee).meas_pattern = meas_pattern(stim_pattern(1:end-1) == 0,:);
end        
        
m = setup_ref_electrodes(m,ref); % Set up reference electrode

%% Set conductivities and run model
printf('Running field simulation ... \n')
if isdeployed, disp('Progress: 25%'), end

if any(named('-sigma')), % Import volume conductivities from fitted transimpedance data
     em = build_forward_model(m,get_('-sigma')); % use custom anisotropic conductivities
else em = build_forward_model(m); % set up non-homogenous anisotropic conductivity
end 

em.fwd_solve.get_all_meas = 1; % Get ALL points
v_meas = fwd_solve(em); % solve EIDORS electric field simulation 
v_meas = v_meas.volt;

if do_PLOT

  ok = cat(1,m.object_id{strncmp(m.object_name,'P_Fa',3)});
  if ~any(ok), 
    ok = cat(1,m.object_id{strncmp(m.object_name,'Fascicle',3)});
  end

  if any(ok)
    ok = unique(m.elems(ok,:));
    scatter3(m.nodes(ok,1), m.nodes(ok,2), m.nodes(ok,3), 10, ...
                      v_meas(ok,1),'o'), colormap(tools.magma)
    caxis(quantile(v_meas(ok,1),[.01 .99]));
  end
end

%% Build output file 

if isdeployed, disp('Progress: 90%'), end


if any(named('-out')), output_name = get_('-out');
else output_name = tools.file('get','out~/extracellular-potential (%d).mat','next');
end

m.name = 'Extracellular potential simulation';
model = m;
v_extracellular = v_meas;

fileparts_output_name = fileparts(output_name);
if ~exist(fileparts_output_name,'dir')
    warning('ViNERS:mesh:makedir','making %s', fileparts_output_name)
    mkdir(fileparts_output_name)
end

printf('Saving %s\n', output_name)
save(output_name,'model','v_extracellular')
         
if isdeployed, disp('Progress: 100%'), end

return 

%% Add one or two reference electrodes to outer surfaces 
function m = setup_ref_electrodes(m,ref) % IN-CONTEXT function

named = evalin('caller','named');
get_ = evalin('caller','get_');
nE = evalin('caller','nE');

if any(named('-reference')),  ref_electrode = get_('-reference');
elseif any(named('-common')), ref_electrode = [2 -1 2]; 
else ref_electrode = [2 -1 2; 1 -1 3]; 
end

d = ref_electrode(1,3); % which side to put on? [xyz] = 1-3
idx = unique(m.boundary(m.boundary_numbers == 1,:));
if ref_electrode(1,2) < 0
     idx(m.nodes(idx,d) > min(m.nodes(idx,d))+10*eps) = []; 
else idx(m.nodes(idx,d) < max(m.nodes(idx,d))-10*eps) = []; 
end

d = ~ismember(1:3,d); 
e_dist = sqrt(sum(m.nodes(idx,d).^2,2)); 
if ~any(e_dist > ref_electrode(1,1))
  warning('pnModel:pointRefElec','COUNTER (ref 1) is a point electrode')
  [~,idx] = min(e_dist);
else idx(e_dist > ref_electrode(1,1)) = []; 
end

m.electrode(ref).nodes = idx;
m.electrode(ref).z_contact = 1; % changing to 0 results in NaN output
m.gnd_node = m.electrode(ref).nodes(1); 

if size(ref_electrode,1) > 1 
  %% Use a second virtual distant electrode for recording reference
  % The mesh-free simplified model suggusted that the "pedestal"
  % which sensitivity peak was sitting on was be artifactual. I
  % confirmed that using a dual-referece model which I am now
  % making the default. 

  d = ref_electrode(2,3); % which side to put on? [xyz] = 1-3
  idx = unique(m.boundary(m.boundary_numbers == 1,:));
  if ref_electrode(2,2) < 0
       idx(m.nodes(idx,d) > min(m.nodes(idx,d))+10*eps) = []; 
  else idx(m.nodes(idx,d) < max(m.nodes(idx,d))-10*eps) = []; 
  end

  d = ~ismember(1:3,d);
  e_dist = sqrt(sum(m.nodes(idx,d).^2,2)); 
  if ~any(e_dist > ref_electrode(2,1))
    warning('pnModel:pointRefElec','GROUND (ref 2) is a point electrode')
    [~,idx] = min(e_dist);
  else idx(e_dist > ref_electrode(2,1)) = []; 
  end

  m.electrode(ref+1).nodes = idx;
  m.electrode(ref+1).z_contact = 1; 
  m.gnd_node = m.electrode(ref+1).nodes(1);

  v_ = @(x) reshape(x,1,[]);
  
  % Update meas pattern
  meas_pattern = sparse(v_((1:nE)' * [1 1]), ... 
                           [1:nE (ref+1)*ones(1,nE)], ...
                            [ones(1,nE) -ones(1,nE)]);
  [m.stimulation.meas_pattern] = deal(meas_pattern); 

  for ii = 1:numel(m.stimulation) % resize sparse stim array
     m.stimulation(ii).stim_pattern(ref+1) = 0;
  end

end

return

%% use EIDORS.mk_image and fill in sigma values from the literature
function em = build_forward_model(fmdl,sigma)

if ~exist('sigma','var')
    sigma(1).name = 'Interstitial';
    sigma(2).name = 'Prostate';      % 0.4243 S/m % http://dx.doi.org/10.1109/TBME.2007.897331
    sigma(3).name = 'PDMS';          % [materials properties] 
    sigma(4).name = 'Fascicle';      % 
    sigma(5).name = 'P_Fascicle';    % https://pubmed.ncbi.nlm.nih.gov/30507555/
    sigma(6).name = 'Epineurium';
    sigma(7).name = 'Blood';         % https://europepmc.org/article/pmc/pmc2914348

    sigma(1).value = 0.66;           % units of S/m
    sigma(2).value = 0.4243;
    sigma(3).value = 1e-12;
    sigma(4).value = [0.570, 0.088, 0.088]; 
    sigma(5).value = 8.7e-4;
    sigma(6).value = 0.083;          % [6.3] ? 
    sigma(7).value = 1.09; 

% Perineurium 0.0021
% Endoneurium
%  longitudinal 0.57
%  transverse   0.083
% Epineurium 0.083
% Saline 2.0
%
% from https://pubmed.ncbi.nlm.nih.gov/14278100/
% THE SPECIFIC IMPEDANCE OF THE DORSAL COLUMNS OF CAT
%
% from http://dx.doi.org/10.1016/0014-4886(65)90126-3
% Specific impedance of cerebral white matter. Nicholson (1965)
%     
%        Normal to fibers    Parallel to fibers
%        Resistive Reactive  Resistive Reactive
% 20 Hz  850 � 150  67 � 32   89 � 40  7 � 8
% 200 Hz 770 � 160  39 � 18   89 � 39  4 � 5
% 2 kHz  770 � 140  35 � 27   78 � 33  6 � 6
% 20 kHz 750 � 150  140 � 110 78 � 33  15 � 22
%
% Table 1: Conductivities used for modeling [12]: 
% AQ Choi, JK Cavanaugh, and DM Durand, "Selectivity of multiple-contact
% nerve cuff electrodes" IEEE Trans Biomed Eng (2001). 
    
end
if ~exist('fmdl','var'), em = sigma; return, end

em = mk_image(fmdl, 0.66);

em.elem_data(:,1,1:3,1:3) = 0;

for ii = 1:length(sigma)
    
  sel = strncmpi(em.fwd_model.object_name,sigma(ii).name, ...
                                 length(sigma(ii).name));
  sel = cat(1,em.fwd_model.object_id{sel});
  
  if ~any(sel), continue, end
  
  if numel(sigma(ii).value) == 3
    em.elem_data(sel,1,1,1) = sigma(ii).value(1);
    em.elem_data(sel,1,2,2) = sigma(ii).value(2);
    em.elem_data(sel,1,3,3) = sigma(ii).value(3);
  else
    em.elem_data(sel,1,1,1) = sigma(ii).value(1);
    em.elem_data(sel,1,2,2) = sigma(ii).value(1);
    em.elem_data(sel,1,3,3) = sigma(ii).value(1);
  end
end

% Check that every element has an assigned conductivity
done = any(em.elem_data(:,:,:),3);
if ~all(done)
    warning('build_forward_model:unsetElements','%d of %d tets have undefined conductivities', ...
        sum(~done), length(done))
end

