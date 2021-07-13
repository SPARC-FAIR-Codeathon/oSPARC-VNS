
function export_H5XP_file(data,varargin)

if nargin == 0, error('Not enough input arguments'), end
if isempty(which('tools.make_SPARC_structure')), path(path,expdir('PN')), end

example = tools.file('get','\\research-cifs.unimelb.edu.au\9560-bi_data\data-bionics\primary\12_005\*.h5xp');
root = h5info(example);

output = 'simulation.h5xp';

if exist(output,'file'), delete(output), end

deep_copy(example,root.Groups(1),output)
deep_copy(example,root.Groups(2),output)

export_EP_data(example,root.Groups(3),data,output)

deep_copy(example,root.Groups(4),output)
deep_copy(example,root.Groups(5),output)


function deep_copy(input,folder,output)


for ii = 1:numel(folder.Datasets)

%     h5create(output_file,ds,size(y),'Datatype','single', ...
%                                   'ChunkSize',[1e3 1],'Deflate',9)
%     h5write(output_file,ds,y * 1e3)
%     h5writeatt(output_file,ds,'StartTime',time_values(1));

  ds = [folder.Name '/' folder.Datasets(ii).Name];
  
  try
    dat = h5read(input, ds );  
  catch C, 
    warning('deep_copy:read_error','%s\n%s: %s', input, ds, C.message), 
    continue
  end
    
  
  if ismember(class(dat), {'double','uint64','uint32','uint16','uint8', ...
                           'single','int64', 'int32', 'int16', 'int8'})

    datatype = {size(dat),'Datatype',class(dat)};     
    h5create(output, ds,  datatype{:} )
    h5write(output, ds, dat)
  else
    % datatype = {1,'TextEncoding','UTF-8'};    
    try
      h5writestr(output, ds, dat)
    catch C
      warning('deep_copy:strwrite_error','%s\n%s: %s', output, ds, C.message), 
      continue

    end
  end
  
  for aa = 1:numel(folder.Datasets(ii).Attributes)
    h5writeatt(output,ds, ...
               folder.Datasets(ii).Attributes(aa).Name, ...
               folder.Datasets(ii).Attributes(aa).Value );
  end
end

for ii = 1:numel(folder.Groups), 
  deep_copy(input, folder.Groups(ii), output)
end


function export_EP_data(input,sample,data,output)

export_root = '/Packed Data/EPs/';


sample = sample.Groups(1).Groups.Groups;


for i0 = 1:numel(data)
  
  trial_name = tools.file('T',data(i0).info(1).path);
  trial_name = sprintf('%s_%s_%s', ...
                        regexp(trial_name,'sub-\d+','match','once'), ...
                        regexp(trial_name,'sam-\d+','match','once'), ...
                        'A'+i0-1);

  pair = data.info.wave_settings(end).bipolar_pair(2,:);
  chan_ = @(s) sprintf('E%dR%d /%s',pair,s);

  var_ = @(n) [export_root trial_name '/' chan_(n)];

  CL = round(data(i0).info.stimulus.current);  
  T_FAKE = CL(round(end*0.6));
  
  h5writestr(output,var_('s_Levels'),sprintf('%d;',CL))    
  
  datatype = {[1 1],'Datatype','single'};
  h5create(output, var_('v_Threshold'),  datatype{:} )
  h5write(output, var_('v_Threshold'), single(T_FAKE) )

  datatype = {size(CL),'Datatype','single'};
  h5create(output, var_('w_Level'),  datatype{:} )
  h5write(output, var_('w_Level'), single(CL) )

  datatype = {size(CL),'Datatype','single'};
  h5create(output, var_('w_Response'),  datatype{:} )
  h5write(output, var_('w_Response'), single(0*CL) )
  
  for ii = 1:numel(sample.Datasets)
    for aa = 1:numel(sample.Datasets(ii).Attributes)
      try
      h5writeatt(output,var_(sample.Datasets(ii).Name), ...
                 sample.Datasets(ii).Attributes(aa).Name, ...
                 sample.Datasets(ii).Attributes(aa).Value );
      catch C, 
        warning('deep_copy:write_data_attr','%s\n%s: %s', input, var_(sample.Datasets(ii).Name), C.message), 
      end
    end
  end
      
  
  t = data.info.time;
  t0 = (data.info.time) >= 0; 
  upsample_ = @(w) interp1(t',w',min(t):(1/200):max(t),'cubic'); 
  
  % for cc = record pairs
  
  pair = data.info.wave_settings(end).bipolar_pair(2,:);
  chan_ = @(s) sprintf('E%dR%d /%s',pair,s);
  
  for ss = 1:size(data.wave,2)
        
    w = [data.wave(end,ss).wave(t0) data.wave(end,ss).wave(~t0)];
    w = upsample_(w);
    
    run_id = strrep(chan_('~001'),'/~',' ');
    
    var_ = @(n) [export_root trial_name '/' chan_(run_id) '/' n];  
    wave_name = sprintf('w_E%dR%d %dCL',pair,CL(ss));
    
    datatype = {size(w),'Datatype','single'};

    h5create(output, var_(wave_name), datatype{:} )
    h5write(output, var_(wave_name), single(w) )

    for aa = 1:numel(sample.Groups(1).Datasets(2).Attributes)
      try
      h5writeatt(output,var_(wave_name{ii}), ...
                 sample.Groups(1).Datasets(2).Attributes(aa).Name, ...
                 sample.Groups(1).Datasets(2).Attributes(aa).Value );
      catch C, 
        warning('deep_copy:write_data_attr','%s\n%s: %s', input, var_(sample.Datasets(2).Name), C.message), 
      end
    end    
  end
end







function h5writestr(filename, dataset, str)
% Matlab's h5write(filename, dataset, data) function doesn't work for
% strings. It hasn't worked with strings for years. The forum post that
% comes up first in Google about it is from 2009. Yeah. This is terrible,
% and evidently it's not getting fixed. So, low level functions. Fun fun.
%
% What I've done here is adapt examples, one from the hdf group's website
% https://support.hdfgroup.org/HDF5/examples/api18-m.html called
% "Read / Write String Datatype (Dataset)", the other by Jason Kaeding.
%
% I added functionality to check whether the file exists and either create
% it anew or open it accordingly. I wanted to be able to likewise check the
% existence of a dataset, but it looks like this functionality doesn't exist
% in the API, so I'm doing a try-catch to achieve the same end. Note that it
% appears you can't just create a dataset or group deep in a heirarchy: You
% have to create each level. Since I wanted to accept dataset names in the
% same format as h5read(), in the event the dataset doesn't exist, I loop
% over the parts of the dataset's path and try to create all levels. If they
% already exist, then this action throws errors too; hence a second
% try-catch.
%
% I've made it more advanced than h5create()/h5write() in that it all
% happens in one call and can accept data inputs of variable size. I take
% care of updating the dataset's extent to accomodate changing data array
% sizes. This is important for applications like adding a new timestamp
% every time the file is modified.
%
%@author Pavel Komarov pavel@gatech.edu 941-545-7573

%"The class of input data must be cellstring instead of char when the
%HDF5 class is VARIABLE LENGTH H5T_STRING.", but also I don't want to
%force the user to put braces around single strings, so this.
if ischar(str), str = {str}; end

%check whether the specified .h5 exists and either create or open
%accordingly
if ~exist(filename, 'file')
  file = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
else
  file = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
end

% set variable length string type
vlstr_type = H5T.copy('H5T_C_S1');
H5T.set_size(vlstr_type,'H5T_VARIABLE');

% There is no way to check whether a dataset exists, so just try to
% open it, and if that fails, create it.
try
  dset = H5D.open(file, dataset);
  H5D.set_extent(dset, fliplr(size(str)));
catch
  %create the intermediate groups one at a time because evidently the
  %API's functions aren't smart enough to be able to do this themselves.
  slashes = strfind(dataset, '/');
  for i = 2:length(slashes)
    url = dataset(1:(slashes(i)-1));%pull out the url of the next level
    try
      H5G.create(file, url, 1024);%1024 "specifies the number of
    catch %bytes to reserve for the names that will appear in the group"
    end
  end

  %create a dataspace for cellstr
  H5S_UNLIMITED = H5ML.get_constant_value('H5S_UNLIMITED');
  spacerank = max(1, sum(size(str) > 1));

  if ndims(str) <= 2 && min(size(str)) == 1, 
       h5_dims = length(str);
  else h5_dims = fliplr(size(str));
  end
  dspace = H5S.create_simple(spacerank, h5_dims, ones(1, spacerank)*H5S_UNLIMITED);

  %create a dataset plist for chunking. (A dataset can't be unlimited
  %unless the chunk size is defined.)
  plist = H5P.create('H5P_DATASET_CREATE');
  chunksize = ones(1, spacerank);
  chunksize(1) = 2;
  H5P.set_chunk(plist, chunksize);% 2 strings per chunk
  dset = H5D.create(file, dataset, vlstr_type, dspace, plist);

  %close things
  H5P.close(plist);
  H5S.close(dspace);
end

%write data
H5D.write(dset, vlstr_type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);

%close file & resources
H5T.close(vlstr_type);
H5D.close(dset);
H5F.close(file);


return

