
% read_dat_file imports text and mixed text/binary files which follow any
%  of my internal formats, typically ending in .(data-type).dat - this
%  mostly allows for relatively straightforward Matlab/Python/etc 
%  interoperability & lightweight data exchange
% 
% As of v0.3, the following files are supported: 
%  [text-only] *.titr.dat     *.im_pot.dat  
%  [binary]    *.ve.dat       *.splines.dat
% 
% Other files can be handled with a minumum of complaining, provided they 
%  consist of the following: BEGIN [settings, axons, threshold, i_] { ... 
%
%  Author:         Calvin Eiber
%  e-mail address: calvin.eiber@unimelb.edu.au
%  Version 0.2     24 August 2019
%          0.3     30 April 2020 (updated splines.dat to include metadata)

function dat = read_dat_file(filename,varargin)


named = @(v) strncmpi(v,varargin,length(v)); 
in_ = @(v) varargin{find(named(v),1)+1};

persistent fp
if nargin == 0 || ~exist(filename,'file')
    
    if nargin > 0, fn = regexprep(filename,'^[^\.]+\.','*.');
    else           fn = '*.dat';
    end
    
    if isempty(fp) || ~ischar(fp), fp = '.'; end
    
    if nargin > 0 && strcmpi(filename,'splines.dat') % shortcut        
        filename = 'R:\Lab-Keast-Osborne\Labdata_CDE\Data\Martin\pelvicNerve.splines.dat';
    else 
        [fn,fp] = uigetfile({fn; '*.dat'},[],fp); 
        filename = [fp fn];
    end
end

%% Set up output data 

 dat.filename = filename;
[dat.path,dat.filename] = fileparts(filename); 

if contains(filename,'.titr.dat')    
    dat.axons = {};
    dat.data = [];
elseif contains(filename,'.im_pot.dat')
    dat.electrode = {}; 
    dat.fascicle = []; 
    dat.potential = [];
    dat.xyz = [];
elseif contains(filename,'.splines.dat') % Special case
    if any(named('index')),
         dat.index = in_('index'); 
         [dat.outline,dat.coeffs] = read_splines_file(filename,dat.index);
    else [dat.outline,dat.coeffs,dat.info] = read_splines_file(filename);
    end
    return
elseif contains(filename,'.ve.dat') % Special case
    dat = read_binary_datfile(filename);
    return
else
    warning('read_dat_file:unknown_filetype',...
            'dat-file specification "%s" unrecognised', dat.filename)
end

dat.settings = {};
dat.notes = {};

%% Read in a text-only data-file, line by line 

mode = 0; line = 0;
fi = fopen(filename,'rt'); 

while ~feof(fi)
    
    in = fgetl(fi); 
    line = line + 1; 
    
    if ~ischar(in), continue, end
    if isempty(strtrim(in)), continue, end
    
    if mode == 0 % Outside of a defined block 
      if ~any(in == '{'), dat.notes{end+1} = in; continue, end        
      switch lower(regexprep(in,'\s*{.*',''))
        case 'settings',  mode = 1; config = {}; 
        case 'axons', mode = 2; 
            axon = struct('f',[],'n',[],'x',[],'y',[]);
            if ~isfield(dat,'axons'), dat.axons = {}; end        
        case 'threshold', mode = 3; I_mA = []; 
            if ~isfield(dat,'data'), dat.data = []; end
        case 'i_',        mode = 4; V_mV = [];
            if ~isfield(dat,'fascicle')
                dat.fascicle = [];
                dat.potential = [];
                dat.xyz = []; 
            end
            % There's contextual data on this line 
            in = str2double(regexp(in,'-?([0-9]*\.)?[0-9]+','match'));
            dat.fascicle = [dat.fascicle; in(1)  ];
            dat.xyz      = [dat.xyz;      in(2:4)];

        otherwise warning('read_dat_file:unknown_datablock', ...
                          'Unrecognised block "%s" at line %d', in, line)
          mode = -1; 
      end % Switch mode-select
      
    elseif mode == 1 % Settings
      if any(in == '}'), mode = 0; 
        dat.settings{end+1} = config';            
        if isfield(dat,'electrode') && isempty(dat.electrode)
          sel = strncmpi(config,'Elec',4);
          elec.name = regexp(config(sel)','Elec[^\s]*','match','once');
          elec.xyz = cellfun(@str2double, regexp(config(sel), ...
                             '(?<=\s)-?\d+(\.\d+)?','match'),'unif',0);
          elec.xyz = cat(1,elec.xyz{:});
          dat.electrode = elec; 
        elseif isfield(dat,'electrode')
            warning('read_dat_file:MultipleSettingsBlocks',...
                    'Multiple settings blocks in %s', dat.filename)
        end
      else config{end+1} = in; %#ok<AGROW>
      end      
      
    elseif mode == 2 % Axons
      if any(in == '}'), % End-of-block
        % axons come out in alphabet order from the titration_sensor
        ax_name = arrayfun(@(f,d) sprintf ('F_%d ax_%d', f, d), ...
                                       axon.f,axon.n,'Unif',0); 
        [~,index] = sort(ax_name);
        axon.f = axon.f(index);
        axon.n = axon.n(index);
        axon.x = axon.x(index);
        axon.y = axon.y(index);
        dat.axons = axon; mode = 0;
        continue
      end
      in = str2double(regexp(in,'-?([0-9]*\.)?[0-9]+','match'));
      axon.f = [axon.f; in(1)];
      axon.n = [axon.n; in(2)];
      axon.x = [axon.x; in(3)];
      axon.y = [axon.y; in(4)];          
      
    elseif mode == 3 % Threshold (from titration sensor in S4L)
      if any(in == '}'), mode = 0; 
        dat.data = [dat.data I_mA(:)];
      else
        in = str2double(regexp(in,'-?([0-9]*\.)?[0-9]+','match'));
        I_mA = [I_mA; reshape(in,[],1)]; %#ok<AGROW>
      end
      
    elseif mode == 4 % I_ as in I_mem -> V_electrode
      if any(in == '}'), mode = 0; 
        dat.potential = [dat.potential; V_mV];
      else
        in = str2double(regexp(in,'-?([0-9]*\.)?[0-9]+','match'));
        V_mV = [V_mV; reshape(in,1,[])]; %#ok<AGROW>
      end
    
    elseif mode == -1 % Unrecognised block
        if any(in == '}'), mode = 0; end        
    else error('Unknown mode "%d", line %d', mode, line)
    end
end

if numel(dat.settings) == 1, dat.settings = dat.settings{1}; end

fclose(fi); 

%%
if ~isempty(dat.notes)
    index = str2double(regexp(dat.notes{1},'(?<=#)\d+','match'));
    if ~isempty(index) && ~isnan(index)
      dat.outline = read_splines_file(index);
    end
end











%% Read in a binary "splines.dat" file generated by make_random_fascicles.m
function [outlines,coeffs,metadata] = read_splines_file(filename,index)

% Parse input behaviour
if isnumeric(filename), index = filename; filename = ''; end
if isempty(filename), 
    filename = 'R:\Lab-Keast-Osborne\Labdata_CDE\Data\Martin\pelvicNerve.splines.dat';
end

fid = fopen(filename,'r'); 

if ~exist('index','var'), index = []; end

outlines = {};
coeffs  = {}; 
metadata = ''; 

while ~feof(fid)

  header = fread(fid,5,'int16'); % one header per section in the file

% Header should contain : 1 header_size_bytes [or 0 if multiple data]
%                         2 n_data_rows
%           			  3 loops_per_row = 5
%       				  4 points_per_loop = 11
% 						  5 sizeof_rows_bytes

  if isempty(index) % get all loops in file   
     input_loop = 1:header(2);
     outline = cell(size(input_loop))'; 
     coeff   = cell(size(input_loop))';
  elseif index > header(2)
    input_loop = 1:header(2);
     outline = cell(size(input_loop))'; 
     coeff   = cell(size(input_loop))';
  else input_loop = 1;
     fseek(fid, (index-1) * header(5),0); % jump to data at specified index
  end
  
  for idx = input_loop

    nF = fread(fid,1,'int16');
    coeff{idx} = fread(fid,2*header(3)*header(4),'double');
    coeff{idx} = reshape(coeff{idx},2,[],header(3));

    for ff = 1:nF
        loop = spline(1:header(4), coeff{idx}(:,:,ff));
        outline{idx}(:,:,ff) = ppval(loop,linspace(1,header(4),65))';
    end
  end

  if numel(coeff) == 1, % Un-box splines
    outline = outline{1}; coeff = coeff{1};
  end
  
  if ~isempty(index) % Specific index requested, return just this from section #1
      outlines = outline;
      coeffs = coeff;
      fclose(fid);
      return
  end
  
  outlines = [outlines {outline}]; %#ok<AGROW>
  coeffs = [coeffs {coeff}];    %#ok<AGROW>
  
  % is this the last section in the file? 
  if header(1) ~= 0
    if ~feof(fid), metadata = strtrim(fread(fid,'*char')'); end
    if numel(coeffs) == 1, % Un-box splines
      outlines = outlines{1}; coeffs = coeffs{1};
    end
    
    fclose(fid); 
    return
  end
  
  outline = cell(0); 
  coeff  = cell(0);
end

warning('read_dat_file:spline:eof', ...
        'Hit end-of-file while reading %s', dat.filename)
fclose(fid);



%% Read in mixed binary/text file generated by export_fields.py
function dat = read_binary_datfile(filename)


 dat.filename = filename;
[dat.path,dat.filename] = fileparts(filename); 

if contains(filename,'.ve.dat')    
    dat.axons = {};
    dat.voltage = []; 
    dat.data = [];
else
    warning('read_dat_file:unknown_filetype',...
            'dat-file specification "%s" unrecognised', dat.filename)
end

dat.settings = {};
dat.notes = {};

mode = 0; line = 0;
fi = fopen(filename,'r'); % mix of binary and text

while ~feof(fi)    
    
    p = ftell(fi); 
    in = fgetl(fi); 
    line = line + 1; 
    
    if ~ischar(in), continue, end
    if isempty(strtrim(in)), continue, end
    
    if ~any(in == '{'), dat.notes{end+1} = in; continue, end        
    switch lower(regexprep(in,'\s*{.*',''))
        case 'axons' % BINARY axon data
            axon = struct('f',[],'x',[],'y',[]);
          % header = "AXONS { fas_id X Y " = 19
            fseek(fi,p + 19,'bof'); 

            % fi.write(struct.pack('1h',len(axons)))
            nA = fread(fi,1,'uint16');
            p = ftell(fi);
            
            % for each axon, [fac_id X Y] as uint16 double double            
            axon.f = fread(fi,[nA 1],'uint16',16);
            fseek(fi,p+2,'bof');             
            axon.x = fread(fi,[nA 1],'double',10);
            fseek(fi,p+10,'bof');
            axon.y = fread(fi,[nA 1],'double',10);
            fseek(fi,-10,'cof'); 
            
            in = fgetl(fi);
            read_err = ~isempty(in); 
            in = fgetl(fi); line = line + 1; 
            read_err = read_err | ~strcmpi(in,'} ENDSEC AXONS');
            
            if read_err
                warning('read_dat_file:bad_binary_axons', ...
                        'Possibly corrupted axon data, line %d', line-1)
            end            
            dat.axons = axon; 

        case 'settings', config = {}; % TEXT settings data
            while ~feof(fi) 
              in = fgetl(fi); 
              line = line + 1;             
              if any(in == '}')
                dat.settings{end+1} = config';
                break
              else config{end+1} = in; %#ok<AGROW>
              end
            end
            
        case 'outline', % TEXT outline data            
            if isfield(dat,'outline')
                 config = dat.outline;
            else config = {};
            end
            
            while ~feof(fi) 
              in = fgetl(fi); 
              line = line + 1;             
              if any(in == '}')
                dat.outline = config;
                break
              elseif isempty(in), continue
              else config{end+1} = reshape( str2double( regexp( ...                  
                                     in, '-?(\d*\.)?\d+','match')), 2,[])'; %#ok<AGROW>
              end
            end
            
        case 've', % BINARY voltage data            
            axon = struct('z',[],'v',[]);            
          % header = "VE { Z Ve " = 10
          
            fseek(fi,p + 10,'bof');
            nZ = fread(fi,1,'uint16'); % # of points in axon
            axon.z = fread(fi,[nZ 1],'float'); % z coordinates
            axon.v = fread(fi,[nZ 1],'float'); % voltages
            
            in = fgetl(fi);
            read_err = ~isempty(in); 
            in = fgetl(fi); line = line + 1; 
            read_err = read_err | in(1)~='}';
            if read_err
                warning('read_dat_file:bad_binary_axons', ...
                        'Possibly corrupted voltage data, line %d', line-1)
            end
            
            if isempty(dat.data), dat.data = axon;
            else                  dat.data(end+1) = axon; 
            end
        otherwise warning('read_dat_file:unknown_binary_datablock', ...
                         ['At line %d: unrecognised block "%s" '    ... 
                          '(possibly containing binary data)'], line, in)
    end % Switch mode-select
end

if numel(dat.settings) == 1, dat.settings = dat.settings{1}; end
fclose(fi); 



if contains(filename,'.ve.dat') % File-type-specific post-processing    
  if length(unique(cellfun(@length,{dat.data.z}))) == 1 % All data is the same size  
       dat.axons.z = dat.data(1).z;
       dat.voltage = [dat.data.v]';
       dat = rmfield(dat,'data');
  else dat = rmfield(dat,'voltage');
  end
end
      
      
      
      
      
      






    