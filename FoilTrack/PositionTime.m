% PositonTime
%   Calculate the position and displacement rate of the foils using displacement data produced by ImageAnalysis script. 
%   It is designed to be used on data which is too noisy to do on a per image
%   basis or for data which has very small smplitudes. 
%
%   syntax: PositonTime(*position_file_name.mat, options)
%     options: 'window', value  -- number of images to use in calculating the strain-rate
%                                       The default value is 5
%              'plot'           -- plot the validation figures.
%              'surface coefficients', m, n -- set the order of the surface
%                                       coefficients for the fitting 
%                                       The default values are 6,1
%
%   See Also: ImageAnalysis

%   $ version: 0.1 $ Feb 2022 $
%       Simon Hunt. 2022

% This version is to make sure the correct length of samples is reported in
% the output file. 


% = PhasesSSD( %<-disp analysis row->%
% = PhasesSSD( %<-disp analysis row->%
% = PhasesSSD( %<-disp analysis row->%

function PositionTime(file_name, varargin)

tic
fprintf(1, '\n Rates in %s \n', file_name);

warning off all

%% get settings and variables.

%set default options
window      = 6;
begin       = -1;
cease       = -1;
plots_on    = 0;
output_save = 1;
overwrite  = 0;
degree = 1;

% prase varargin...
iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'window'
            window = varargin{iarg+1};
            iarg = iarg + 2;
        case 'degree'
            degree = varargin{iarg+1};
            iarg = iarg + 2;
        case 'begin'
            begin = varargin{iarg+1};
            iarg = iarg + 2;
        case 'cease'
            cease = varargin{iarg+1};
            iarg = iarg + 2;
        case 'plot'
            plots_on = 1;
            iarg = iarg + 1;
        case 'overwrite'
            overwrite = 1;
            iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end


% get experiment options and if not then make them.
expt = FileTitleInformation(file_name);
varfilename = char(strcat(expt(1), '.mat'));

%list of variables need to get/generate.
filevariables = {'length_type', 'get_rid', 'symm_type', 'scaling', 'search', 'analysis_type'};

%reads the experiment analysis options file. -- the file is called 'expt_name'.mat
if exist(varfilename, 'file')
    vars = load(varfilename, filevariables{:});
end
if exist("search","var")
    vars.search=search;
end

% if the variables are not read from 'expt_name'.mat run AnalysisOptions to generate them.
if sum(isfield(vars,filevariables)) ~= (length(filevariables))
    [out] = AnalysisOptions;
    load(varfilename, filevariables{:});
end

%make filename for output file
root = file_name(1:end-19);
%root_length = strfind(root,'_SSD');
%root = root(1:root_length-1);
tail = '_displacements_rates';
%out_file_name1 = [root,'_window_positions.txt']; %removes 'position_change' from the end of the input data file name.
%out_file_name2 = [root,'_window_lengths.txt']; %removes 'position_change' from the end of the input data file name.
out_file_name3 = [root,'window_PositionChangeRate.txt']; %removes 'position_change' from the end of the input data file name.
out_file_name4 = [root,'window_LengthChangeRate.txt']; %removes 'position_change' from the end of the input data file name.


%% data setup 

%get box reference positions
if strcmpi(file_name(end-18:end),'position_change.txt')
    boxes = load([file_name(1:end-19), 'boxes']);
else
    boxes = load([file_name, '_boxes']);
end

%import the postion change data
[foil_change_data_x, image_time, box_positions_x, number_boxes_x, ...
    number_images_x] = ReadPositionChangeFile([file_name]);
%fix reference lengths from rounded values in file
box_positions_x=[boxes.boxY';boxes.boxX'];


% Smooth the time stamps (if necessary) and make times start from 0
arg_to_pass = {}; if plots_on == 1, arg_to_pass = {arg_to_pass{:}, 'plot'}; end
image_time = TimeStampSmooth(image_time, arg_to_pass{:}); %required for tiff image based data sets. If netcdf based the returns same array as input.
%image_time = repmat(image_time',dat_size(1),1); %makes image_time array equivalent size to other arrays


%debugging function.
if (0)
    fprintf('Size of data array      : %d %d %d \n', size(pos_array,1), size(pos_array,2), size(pos_array,3))
    fprintf('Size of time stamp array: %d %d    \n', size(image_time,1), size(image_time,2))
    keyboard;return
end

%% positions
if vars.analysis_type == 'ref'
    position = foil_change_data_x + mean(boxes.boxY');
else
    error('this option is not implemented')
end

%% rate shape for individual boxes
disp('Fit polymonial shape for each foil')

if begin == -1
    process_start = 1;
else 
    process_start = begin;
end
if cease == -1 || cease > length(image_time)
    processs_end = length(image_time);
else 
    processs_end = cease;
end
if window == -1
    window = processs_end-process_start+1;
end

for j = process_start : processs_end-window+1

    pos_now  = foil_change_data_x(j:j+window-1,:);
    time_now = image_time(j:j+window-1)-min(image_time);
    
    for i = 1: number_boxes_x
    
        pp = polyfit(time_now,pos_now(:,i),degree);

        pder = polyder(pp);



        %solve for position at average time
        mean_time_now = mean(time_now);

        box_rate(j,i) = polyval(pder, mean(time_now));


        if 0%plots_on==1

            subplot(2,1,1)
            plot(min(X):(max(X)-min(X))/200:max(X), multi_displacements, '-')
            hold on
            plot(unique(X), polyvalnm_solve2(surf_diff, 0, unique(X)),'o')
            xlabel('Relative time (s)')
            ylabel('SSD Displacement surface minimum')

            subplot(2,1,2)
            plot(min(X):(max(X)-min(X))/200:max(X), polyval(slope,min(X):(max(X)-min(X))/200:max(X)), '-')
            hold on
            plot(unique(X), polyval(slope,unique(X)),'o')
            xlabel('Relative time (s)')
            ylabel('Gradient of surface minimum')
            
            %keyboard
        end

    end
    refs{j} = NaN;
    image_ids{j} = j;
    timestamp(j) = mean_time_now;%t_mean(j);
        
    
    %if required plot data
    if 0%plots_on == 1
        %surface_plot2(X, pos_now, SSD_now, squeeze(fits(j,:,:)), NaN)
        pause
    end
end

%figure
%^plot(slope_at_mean)

%% combine data and prepare for output
l  = position(:,2:end) - position(:,1:end-1);
dldt = box_rate(:,2:end) - box_rate(:,1:end-1);


% figure
% subplot(2,1,1)
% plot(l)
% subplot(2,1,2)
% plot(dldt)



%% generate output header contents

% load experiment options.
try
    expt = FileTitleInformation(file_name);
    varfilename = fullfile(cd,[expt{1},'.mat']); 
catch %#ok<CTCH>
    expt = '';
    varfilename = file_name(1:6);
%     [~,varfilename,~] = fileparts(file_name);
end
if isempty(expt) == 1
   [~, expt,~] = fileparts(file_name);
    varfilename = [expt,'.mat'];
end
%list of variables need to get/generate.
filevariables = {'min_type', 'search', 'analysis_type', 'spot_removal', 'NaN_bg', 'out_data_type', 'expt_location', 'offset'};
%reads the experiment analysis options file. -- the file is called 'expt_name'.mat
if exist(varfilename, 'file') == 2
    AnalysisVariables = load(varfilename, filevariables{:});
else %option was added for images not conforming to X17B2 nomenculture (i.e. DESY experiments)
    % tries to load filevariables from all *.mat files. If it finds them then it is fine.
    
    mat_files = dir('*.mat');
    for x = 1:size(mat_files,1)
        AnalysisVariables = load(fullfile(cd,mat_files(x).name), '-mat', filevariables{:});
    end
end
% if the vaibles are not all there makes them as blanks and passes them to AnalysisOptions to fill in.
if ~exist('AnalysisVariables', 'var')
    AnalysisOptions(varfilename);
    AnalysisVariables = load(varfilename, filevariables{:});
    AnalysisVariables.analysis_type;
end
AnalysisVariables.begin = begin;
AnalysisVariables.cease = cease;
AnalysisVariables.window = window;
%AnalysisVariables.surface_coef = [m,n];

%get the box positions
Boxes = boxes
if isfield(Boxes, 'Image_rotation') == 1
    Rotation = Boxes.Image_rotation;
else
    Rotation = 0;
end

%create array to write in both positionchange headers and SSD files
headers.file_name = file_name;
file_information = FileTitleInformation(file_name);
headers.run_name  = file_information{1};
headers.caller    = [mfilename('fullpath'),'.m'];
headers.variables = AnalysisVariables;
headers.boxes     = Boxes;

%get the reference ID. 
headers.refid     = -1;%data.refid;
        

%% write output files
for x = 3:4

    %select the file and the data
    if x==1
        %position changes per foil
        of = out_file_name1;
        v = position;

    elseif x==2
        %length changes
        of = out_file_name2;
        v = l;

    elseif x==3
        %position change rate
        of = out_file_name3;
        v = box_rate;

    elseif x==4
        % length change rate
        of = out_file_name4;
        v = dldt;

    else
        disp("oops")
    end

    %do the file writing
    if exist(of,'file') == 0 || overwrite==1
        status = WritePositionChange(of, 'OpenNew');
    else
        status = WritePositionChange(of, 'Open');
    end
    status = WritePositionChange(of, 'header', headers);

    %write the positions
    for j = process_start : processs_end-window+1
        %write the offsets to the output file.
        status = WritePositionChange(of, 'values', v(j,:), refs(j), image_ids(j), timestamp(j)+min(image_time(:)));
    end
    WritePositionChange(of, 'close');

end

