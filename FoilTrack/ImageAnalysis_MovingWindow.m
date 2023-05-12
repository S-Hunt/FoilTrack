%% ImageAnalysis_Loop.
%   Loops through experiment data files to use each image as a reference for the displacements.
%   Ideally will be able to pair each image once but for now uses a brute force approach.
%
%   Designed and used on Zn_06. Not tested on other experiments. 
%   
%   Syntax:
%       ImageAnalysis_Loop_devCu02('Cu_02_20tons_ramping_All.lst', 'boxes',...
%               'exclude', 'Cu02_all_Exclusions.txt', ...
%               'window', 60, ...
%               'step', 60, ...
%               'IA Options', {'Dark Field', 'Cu02_DarkFeildImage.tif'})
%

% Simon Hunt 
% May 2016 - 2018

% images recorded every 1.5 seconds
% apparent period from data 30s.
% therefore we have 20 frames / cycle. 

function ImageAnalysis_RefLoop_dev(varargin)


% approximately 20 frames / cycle
global data_points xyz  

%FIX ME: Need a switch be be able to run the analysis forwards or
%backwards. over the xyz loop

%process varargin and defaults
data_points = 81;
direction = 'forward';
ImFunctToCall = str2func('ImageAnalysis');
step = round((data_points-1)/2);

Func_opts = {};
Excl_list = [];
pwr = 0;
%process varargin
file_str = varargin{1}; %files to work on e.g. '*.nc'
IA_funct = varargin{2}; %disp or boxes or alternatives
iarg = 3;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'forward', 'backward'}
            direction = varargin{iarg};
            iarg = iarg + 1;
        case {'window'}
            data_points = varargin{iarg+1};
            iarg = iarg + 2;
        case {'step'}
            step = varargin{iarg+1};
            iarg = iarg + 2;
        case {'function'}
            ImFunctToCall = varargin{iarg+1};
            iarg = iarg + 2;
        case 'exclude'
            Excl_list = varargin{iarg+1};
            iarg = iarg + 2;
        case 'power' %list powers rather than temperature in names.
            pwr = 1;
            iarg = iarg + 1;
        case 'ia options'
            Func_opts = varargin{iarg+1};
            iarg = iarg + 2;
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end


%Get list of files to be processed. (should be *.nc or *.lst files)
files = dir(file_str);
% file_names = {files.name};


%% generate a position change / SSD file for each reference image.

%get the number of images in the file.
times = ImFunctToCall(file_str, 'contents', 'timestamps');
array_size = length(times); %proxy for the number of images in the file.

r_name = ImFunctToCall(file_str, 'contents', 'run_name');

%start position
if strcmpi(direction, 'backward') == 1
    st = array_size-data_points+1;
    en = 1;
    step = step * -1;
else
    st = 40; %1 %%FIX ME This is for the Zn06 window loop. Needs fixing later
    en = array_size-data_points+1;
    step = step;
end

%process exclusions if they exist.
if isnumeric(Excl_list) || isempty(Excl_list)
    %exclusions are OK. Cary on. nothing else to go here
    exclud = Excl_list;
else
    try %atempt to open a file. This should be the only other option.
        excl = importdata(Excl_list);
        
        if numel(excl) == 1
            %the numbers in the exclusion file are relative to the start of
            %the data series not the numbers stamped in the file names.
            exclud =  eval(['[',excl{1},']']);
        elseif numel(excl) == 2
            %the numbers in the exclusion file are the numbers stamped in the file names.
            exclud =  eval(['[',excl{2},']']);
            
            %the list start number is in the first row of the file after
            %'start'
            lst_start = str2num(excl{1}(strfind(excl{1}, ' ')+1:end));
            
            %adjust the exclusion list to match the positions in the list.
            exclud = exclud-(lst_start-1);
            
        end
    catch
        %this just passes the name of the input file to the exclusions
        %function -- used for the Zn_02 experiment.
        
        exclud = files;
    end
end

%loop over the file
for xyz = st:step:en %array_size-data_points+1:-1:st
    
    %check exclusions.
    exclude = exclusions_list(exclud, xyz, data_points);
    
    if exclude == 0
        
        % set new variables
        s = xyz;
        e = xyz + data_points - 1;
        
        %new reference image
        ref_im_num = xyz + ((data_points-1)/2);
        
        run_name =  ([r_name,'_ref', sprintf('%03d',ref_im_num)]);
        
        %change the temperature in the file name
        middle_time = mean(times(s:e));%mean(times);
        temperature = temp_time(middle_time, times(s), times(e), file_str, pwr); %return temperature as a string inless pwr == 1 then power.
        if ~isnan(temperature)
            %temp_str = sprintf('%04.2fC', temperature);
            %temp_str = strrep(temp_str, '.', 'pt');
            run_name = strrep(run_name, 'ramping',temperature);
        end
                
        %set the switch for fast boxes
        if xyz == st% first point to be processed
            %...array_size-data_points+1;%e == data_points % last point in file -- first to be processed.
            fast = 0;
        else
            fast = 1;
        end
        
        %Run the function. ImageAnalysis or one of its relatives.
        ImFunctToCall(file_str, IA_funct,...
            'ref_id', ref_im_num,...
            'start', xyz,...
            'max_frames', e,...
            'run_name', run_name,...
            'fast_box', fast,...
            Func_opts{:});%,...
        
    elseif exclude == 1
        disp(['Excluding ',files.name, ', image ', num2str(xyz), ' as a reference image.'])
    else
        error('Damn this has gone wrong')
    end
end


end

function out = exclusions_list(Excl_list, step, step_range)

%this is the list of images in the Zn06 experiment which have no data in
%them -- i.e. diffraction was being done or the beam was down.
%listed as frame numbers in the data files. 
% persistent exclude_list


if ~isnumeric(Excl_list) && ~isempty(Excl_list)
    %the exclusion lists should be numeric or empty. Zn02 is the historic
    %exception. Its exclusions are listed here -- but should be put into
    %exclusion files.
   
    switch Excl_list.name
        case {'Zn_06_27tons_25C_30s_001.nc'
                'Zn_06_27tons_NoSinePump_150C_002.nc'
                'Zn_06_27tons_ramping_30s_003.nc'
                'Zn_06_27tons_ramping_30s_007.nc'
                'Zn_06_27tons_ramping_30s_008.nc'
                'Zn_06_27tons_ramping_30s_009.nc'}
            Excl_list = [];
            
        case 'Zn_06_27tons_ramping_30s_004.nc'
            Excl_list = [159:247, 446:555];
            
            
        case 'Zn_06_27tons_ramping_30s_005.nc'
            Excl_list = [446:560, 612:727,1393:1587,1867:1963,2132:5000];
            %5000 is well beyond the end of the file and is just a filler
            %number
            
            
        case 'Zn_06_27tons_ramping_30s_006.nc'
            Excl_list = [1:93, 907:1093];
            
        otherwise
            
            warning('Filename not recognised')
            Excl_list = [];
    end
end

step_range = step:step+step_range;
if sum(ismember(step_range, Excl_list))>=1
    out = 1; %exclude
else
    out = 0; %keep
end
        
end




function temp_out = temp_time(time, start_t, end_t, lst_file, pwr)

%deermine file type
[~,~,ext] = fileparts(lst_file);

switch ext
    case '.lst'
        
        %get times from lst file
        t = Image_Functions_Tiff('times', lst_file);
            
        if pwr == 0
            %get temperatures from list file
            T = Image_Functions_Tiff('temperature', lst_file);
            
            s = find(start_t == t);
            e = find(end_t == t);
            
            temp_out = mean(T(s:e));
            
            temp_str = sprintf('%04.2fC', temp_out);
            
        elseif pwr == 1
            %get temperatures from list file
            T = Image_Functions_Tiff('power', lst_file);
            
            s = find(start_t == t);
            e = find(end_t == t);
            
            temp_out = mean(T(s:e));
            
            temp_str = sprintf('%04.2fW', temp_out);
        end
            
        temp_str = strrep(temp_str, '.', 'pt');
        
    case '.nc'

        %Zn_07 experiment
        %time stamps of NETcdf file and tempertures from
        %lab notebook.
        time_temp = [1709422	150
            1709874	163
            1711296	200
            1711546	207
            1711999	219
            1712945	245
            1713202	255
            1713895	277
            1714140	285
            1715104	316
            1715290	323
            1715588	332
            1716009	347
            1716418	362
            1717800	413
            1717894	416
            1718180	425
            1718364	433
            1718899	449
            1719904	483
            1721544	536
            1722033	548
            1722370	558];
        
        mean_time = mean(time_temp(:,1));
        std_time = std(time_temp(:,1));
        poly_time = polyfit((time_temp(:,1)-mean_time)/std_time, time_temp(:,2)-time_temp(1,2) , 3);
        
        if time<time_temp(1,1)
            temp_out = 150;
        else
            temp_out = polyval(poly_time, (time-mean_time)/std_time)+time_temp(1,2);
        end

        temp_str = sprintf('%04.2fC', temp_out);
        
    otherwise
        error('File type not recognised')
end

temp_out = strrep(temp_str, '.', 'pt');

end



