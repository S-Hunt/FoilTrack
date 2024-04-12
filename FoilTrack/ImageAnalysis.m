% ImageAnalysis
%  Thermal diffusivity analysis driver code for Image file (primarily NetCDF) data
% 
% Read images for strain analysis from a NetCDF file or list file, calculate
% missfit and write to file for further analysis.
% 
%   Syntax:
%   ImageAnalysis('FILE_NAME','TO_DO', 'options'...)
%
%   Arguments:
%   1: FILE_NAME - name of file to read. Also used as seed for output file name.
%   2: TO_DO - tells the routine which process to do. The options are:
%         - 'parse'   -- reads all the files so that windows knows where data is for subsequent analysis
%         - 'tiffs'   -- makes tiff images of the ref images in the files
%         - 'boxes'   -- makes boxes for displacement analysis
%         - 'disp'    -- disaplcement analysis
%         - 'view'    -- displays the images sequentially for viewing
%         - 'movie'   -- turn the image data into a .avi movie
%         - 'profile' -- makes a profile through band of all images in experiment and then finds minima of each foil in each frame
%         - 'scale'   -- takes a suite of images in which the press is moved and calculates the scale for the images in microns/pixel. 
%   3: options - this is a list of varaiable names and their new values. It
%   is used or overriding the set options for values in the analysis. Use
%   with Caution.
%
%   See also AnalysisOptions, FileTitleInformation, SelectBoxes, displacement_analysis
%
%   Andrew Walker 2009, 2012 / Simon Hunt 2009-2019
%   $ version: 2.4 $ 21st September 2019 $
%       - for details of changes between versions see end of the file.


function varargout = ImageAnalysis(file_name,to_do,varargin)

tic %start timing for file.

if nargin < 2
    error('There are not enough input arguments. See file help for options.')
end


%% Find output-files location and names
%make run name for outfiles and file type.
[~, run_name, file_type] = fileparts(file_name);

% defines which functional to use based on file type
if strcmpi(file_type,'.nc') == 1
    Im_procs = @Image_Functions_NetCDF;
elseif strcmpi(file_type,'.lst') == 1
    Im_procs = @Image_Functions_Tiff;
elseif strcmpi(file_type,'.edflist') == 1 | strcmpi(file_type,'.edflst') == 1
    Im_procs = @Image_Functions_EDF;
elseif strcmpi(file_type,'.nxslst') == 1
    Im_procs = @Image_Functions_NXSlist;
else
    error('The file type is not recognised')
end

% Where are the files and where do we want the output file to go?
where = cd;


%% load experiment options.
try
    expt = FileTitleInformation(file_name);
    varfilename = fullfile(where,[expt{1},'.mat']); 
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
        AnalysisVariables = load(fullfile(where,mat_files(x).name), '-mat', filevariables{:});
    end
end
% if the vaibles are not all there makes them as blanks and passes them to AnalysisOptions to fill in.
if ~exist('AnalysisVariables', 'var') | length(fieldnames(AnalysisVariables)) ~= (length(filevariables))
    AnalysisOptions(varfilename);
    AnalysisVariables = load(varfilename, filevariables{:});
    AnalysisVariables.analysis_type;
end

%% load images + setup

% opens image file(s)
file_id = Im_procs('open',file_name);

% get timestamps and unique ID's from data file(s)
timestamps = Im_procs('times',file_id);
uids = Im_procs('unique',file_id);

%get number of images and dimensions of the image array
[max_frames image_size(1) image_size(2)] = Im_procs('NumSize',file_id);

%determine reference image number and start position from analysis type.
% CAUTION: Do not change these lines without making corresponding changes
% in PhasesSSD/'combine data and prepare for output'. The position of the
% reference image is needed to correct the lengths for the phase angle of
% the reference.
if strcmpi(AnalysisVariables.analysis_type,'ref')
    start = 1;
    ref_id = round(length(uids)/2);%+15; %FIX ME THIS +15 IS ADDED FOR THE 16IDB EXPERIMENTS!
elseif strcmpi(AnalysisVariables.analysis_type,'seq')
    start = 2;
    ref_id = 1;
else
    error('The requested process type is not recognised')
end

%switch for making lots of boxes very quickly
fast_box = 0; %default is off.

%for plots in displacement analysis
plots_on = 0;
        
%set switches for data types to be saved. 
%assume both will be saved. 
save_SSD = 1;
save_disp = 1;
if strcmpi(AnalysisVariables.out_data_type, 'Displcements') == 1
    % if displcamcents turn off SSD
    save_SSD = 0;
elseif strcmpi(AnalysisVariables.out_data_type, 'SSD') == 1
    % if SSD turn off displcamcents
    save_disp = 0;
end

%default for dark field subtraction -- 0 ie. off
Subtract_DF = 0;

%% override defaults (process varargin)
% This is a dangerous option to override variables set in the code before this point. 
% N.B. cannot override "BlockFrames" using these lines - but it should not
% be need to do so because one can over ride 'max_frames' and 'start' to
% the same effect. 
if strcmpi(to_do, 'contents') ~= 1 %i.e. if we are not looking to output an array. 
    variables = who;
    Analysisfields = fields(AnalysisVariables);
    
    iarg = 1;
    while iarg <= (length(varargin))
        
        if ismember(varargin{iarg},variables)
            %this function may be a better way of forcing the value
            %changes.
%             assignin('caller', varargin{iarg}, varargin{iarg+1});
            
            %change the variable in varargin to its value set by the next entry in varargin
            if isnumeric(varargin{iarg+1})
                eval([varargin{iarg},' = ', sprintf('%20.8f',varargin{iarg+1}),';']);
                
            else
                eval([varargin{iarg},' = ''', varargin{iarg+1},''';'])
            end
            iarg = iarg + 2;
            
        elseif ismember(varargin{iarg},Analysisfields)
            %change AnalysisVariables setting.
            Analysisfields.(varargin{iarg}) = varargin{iarg+1};
            iarg = iarg + 2;
            
        elseif strcmpi(varargin{iarg},'Dark Field')
            
%             DF_image = imread(varargin{iarg+1});
            %assume the dark field image is tif.
            DF_image = Image_Functions_Tiff('getimage',varargin{iarg+1}, AnalysisVariables.expt_location, 1, 1, image_size);
            Subtract_DF = 1;
            
            iarg = iarg + 2;
        elseif strcmpi(varargin{iarg},'FFT Dark Field') | strcmpi(varargin{iarg},'FFTDarkField')
            
%             DF_image = imread(varargin{iarg+1});
            %assume the dark field image is tif.
            DF_image = Im_procs('getimage',varargin{iarg+1}, AnalysisVariables.expt_location, 1, 1, image_size);
            %DF_image = Image_Functions_Tiff('getimage',varargin{iarg+1}, AnalysisVariables.expt_location, 1, 1, image_size);
            DF_image = imbinarize(DF_image,.9);
            if isequal(size(DF_image), image_size) %catch error from FFT which makes the filter. This is due to how the images are read into FFT.
                if isequal(size(DF_image), fliplr(image_size)) %flip DF_image
                    DF_image = DF_image';
                else
                    error('The image filter is not recognised as the right size')
                end
            end
            Subtract_DF = 2;
            
            iarg = iarg + 2;
        else
            error(['Unknown variable or field: ' varargin{iarg}]);
        end
    end
end
   

%outfile names
outfile_name =  fullfile(where,[run_name,'_position_change.txt']); %file which the offsets are written to
box_file_name = fullfile(where,[run_name,'_boxes.mat']); %the file in which the box positions are stored
SSD_file_name = fullfile(where,[run_name,'_SSD.mat']); %the file in which the box positions are stored

%get reference image (called image1)
image1 = Im_procs('getimage',file_id, AnalysisVariables.expt_location, ref_id, 1, image_size);
if Subtract_DF == 1
    image1 = image1-DF_image;
elseif Subtract_DF == 2
    image1 = FFTImageFilter(image1, DF_image);
end
% imagesc(image1), colorbar, keyboard


    
%% process
switch to_do

    case 'parse'
        %there is nothing in here delibarately. All this option is for is to get matlab/windows to parse the file structure.
        %This only needs to be done every time windows is restarted, not when matlab is restarted.
        
        %file_id is closed at the end of the switch/case function
        
    case 'contents'
        % output array from the script. 
        if numel(varargin) ~= 1
            error('Too many outputs requested for ''contents''.')
        end
%         Im_procs('close',file_id);
        varargout = {eval(varargin{1})};
%         return
        
        %file_id is closed at the end of the switch/case function

    case 'export'
        % extracts all the images from a multi-file format. 
        % it will also convert the image sequence into tiffs if called. 
        times_all = [];
        block_frames = 200;
        for j = start : block_frames : max_frames
            
            if j+block_frames > max_frames
                block_frames = max_frames - j +1;
            end
            
            Im = Im_procs('getimage',file_id, AnalysisVariables.expt_location, j, block_frames, image_size,'original');
            im_times = Im_procs('times',file_id, AnalysisVariables.expt_location, j, block_frames, image_size);
            times_all = [times_all; im_times];

            for k = 1:block_frames
                %make_tiffs(Im(:, 450:450+260,k),run_name, k+(j-1)*block_frames);
                make_tiffs(Im(:, :,k),run_name, k+(j-1)*block_frames);

            end
        end
        
        %end time (time the last file in the squence was modified)
        details = dir(file_name);
        t0 = datenum(details.date);
        
        %file_id is closed at the end of the switch/case function
        %t0 = datenum(2020,0,1,0,0,0);
        for i =1 : max_frames
            %set time
            format = 'mm/dd/yy HH:MM:SS.FFF';
            t = (times_all(i) - times_all(1))/24/60/60;
            times_str = datestr(t0+t, format);
            %set filename 
            fnam = make_tiff_names(run_name, i);
            
            %change the creation date of the file using bash function !SetFile 
            % An example of the command is:
            % !SetFile -d "12/31/2000 23:59:59" Zn_08_27tons_117C_100s_027_image00001.tif
            exec_str = sprintf("!SetFile -d ""%s"" %s", times_str, fnam);
            disp(exec_str)
            eval(exec_str)

            %set the modified and accessed times/dates for the extracted images
            % An example of the command is:
            % touch -a -m -d 2022-12-31T11:22:33.5 Zn_08_27tons_117C_100s_027_image00001.tif
            format = 'yyyymmddHHMM.SS';
            times_str = datestr(t0+t, format);
            exec_str = sprintf("!touch -a -m -t %s %s", times_str, fnam);
            disp(exec_str)
            eval(exec_str)
               
            %keyboard
                
        end

       
    case 'reference tiffs'
        make_tiffs(image1,run_name, ref_id);
        
        %file_id is closed at the end of the switch/case function

        
    case 'DarkField'
        FFT(image1, run_name);
        
        return
        
        %file_id is closed at the end of the switch/case function
        
    case 'boxes'
        % checks to see if the box position file exists and if not generates it.
%         SelectBoxes_rewrite(image1, outfile_name, run_name, ref_id, length(uids));
        SelectBoxes(image1, outfile_name, run_name, ref_id, length(uids), fast_box);
%         Im_procs('close',file_id);
%         return %-- this is here so that the box files can be batch generated.

        %file_id is closed at the end of the switch/case function
        
    case {'view', 'profile', 'movie', 'scale', 'screen', 'profile track'}
        global x_im y_im %remeber the positions of the profile selections
        
%         clf
        block_frames = 200;
        frame_step = 1;
        
        %loads the box positions.
        if exist('box_file_name', 'file')
            Boxes = load(box_file_name);
            if isfield(Boxes, 'Image_rotation') == 1
                Rotation = Boxes.Image_rotation;
            else
                Rotation = 0;
            end
        else
            Rotation = 0;
        end
        
        
        
        if strcmpi(to_do, 'movie') == 1
                % Create the output movie
            % Note that on a Mac this needs the (legacy) Quicktime 7 to play,
            % not Quicktime 10.x which is present by default. Apple have
            % dropped the codec used by Matlab...
            % Avoid having to store the whole file in memory (c.f. movie2avi)

            fps = 1/(timestamps(2) - timestamps(1)) *10;
            frame_step = 20;
           
            
            writerObj = VideoWriter([file_name(1:end-3), 'avi']);
            open(writerObj);
%             aviobj = avifile([file_name(1:end-3), '.avi'],'FPS', fps) ;
        elseif strcmpi(to_do, 'view') == 1
            if nargin >= 3
                if strcmpi(varargin{1}, 'start') == 1
                    start = varargin{2};
                else
                    start = 1;
                end
            end
        elseif strcmpi(to_do, 'profile') || strcmpi(to_do, 'profile track')
            start = 1;
        end
        
        hij = figure;
        for j = start : block_frames : max_frames
            
            if j+block_frames > max_frames
                block_frames = max_frames - j +1;
            end
            
            Im = Im_procs('getimage',file_id, AnalysisVariables.expt_location, j, block_frames, image_size);
                    
            %rotate image if required
            if Rotation ~= 0
                Im = rotate_images(Im, Rotation);% Im_procs('rotate',image1, Rotation);
            end
            
            % remove dark field if required.
            if Subtract_DF == 1
                Im = Im-repmat(DF_image,[1,1,size(Im,3)]);
            elseif Subtract_DF == 2
                Im = FFTImageFilter(Im, DF_image);
                
                %[~,diff] = FFTImageFilter(Im(:,:,1), DF_image);
                %Im = Im-repmat(diff,[1,1,size(Im,3)]);
            end
            
            profiles_out=[];
            for x = 1 : frame_step : block_frames
                framenum = j + x - 1;
                image = Im(:,:,x);
                
                limitbot = prctile(image(:), 0);
                limittop = prctile(image(:), 100);
                clims = [limitbot limittop];
                %imagesc(image); axis square; colormap(gray); 
                imagesc(image,clims); axis square; colormap(gray); %line inserted to make vidoe of spt_004 experimental data
                title([num2str(framenum), ' / ', num2str(max_frames)]);
                
                if strcmpi(to_do, 'profile') == 1 || strcmpi(to_do, 'profile track') == 1 ||strcmpi(to_do, 'scale') == 1
                    if exist('y_im','var') == 0 | isempty(y_im)==1
                        [x_im,y_im] = ginput(2);
                        x_im = round(x_im);
                        y_im = round(y_im);
                        
                        set(gcf, 'Visible', 'off')
                        h = waitbar(0,'Calculating profiles...');
                    end
                    prof(framenum,:) = nanmean(image(y_im(1):y_im(2),x_im(1):x_im(2)),2);
                    
                    profiles_out(end+1).profile = prof(framenum,:);
                    profiles_out(end).timestamp = timestamps(framenum);
                    profiles_out(end).uid = uids(framenum);
                    
                    h = waitbar(framenum / max_frames);
                elseif strcmpi(to_do, 'movie') == 1
                    F = getframe(hij);
                    writeVideo(writerObj,F);
%                     aviobj = addframe(aviobj,F);
                elseif strcmpi(to_do, 'screen') == 1
                    
%                     discard = 
                else
                    pause
                end
            end
        end

        if strcmpi(to_do, 'profile') == 1 || strcmpi(to_do, 'profile track') == 1
            
            close(h)
%             figure
            imagesc(prof'), colormap(gray); hold on
            title('Profile through band in images');
            xlabel('image number in series')
            ylabel('mean intensity across profile') 

            if strcmpi(to_do, 'profile track') == 1 
                %select foils/minima to be tracked.
                [x_mins, y_mins] = ginput;
                x_mins = round(x_mins);
                y_mins = round(y_mins);
                for x = 1:length(x_mins)

                    %find local minimum
                    local_y_min = y_mins(x)-50; if local_y_min < 1, local_y_min = 1; end
                    local_y_max = y_mins(x)+50; if local_y_max > size(prof,2), local_y_max = size(prof,2); end
                    [I]=polymin(prof(x_mins(x),local_y_min:local_y_max));
                    I = I + local_y_min;
                    prof_mins(x_mins(x),x) = I;

                    frame = x_mins(x);
                    moving_min = I;
                    window = 20;
                    while frame <= max_frames-1 %finds minima in increasing frames from selected point
                        frame = frame+1;

                        local_y_min = round(moving_min)-window; if local_y_min < 1, local_y_min = 1; end
                        local_y_max = round(moving_min)+window; if local_y_max > size(prof,2), local_y_max = size(prof,2); end

                        [II]=polymin(prof(frame,local_y_min:local_y_max));
                        moving_min = II + moving_min-window;
                        prof_mins(frame,x) = moving_min;

                    end

                    frame = x_mins(x);
                    moving_min = I;
                    while frame >= 2  %finds minima in decreasing frames from selected point
                        frame = frame-1;

                        local_y_min = round(moving_min)-window; if local_y_min < 1, local_y_min = 1; end
                        local_y_max = round(moving_min)+window; if local_y_max > size(prof,2), local_y_max = size(prof,2); end

                        [II]=polymin(prof(frame,local_y_min:local_y_max));
                        moving_min = II + moving_min-window;
                        prof_mins(frame,x) = moving_min;

                    end
                    plot(1:max_frames, prof_mins(:,x),'r.')
                end
            end
            saveas(gcf,[run_name,'_profile_series.jpg'])
            
            
            for x = 1 : length(uids) 
                out_ref{x} = uids(x); %N.B. out_ref and out_nos need to be cell arrays inorder to cope with any change in length of the strings.
            end
            %write output file.
            %fid = fopen(strcat(expt,'_profiles.xls'), 'w');
            %fprintf(repmat('%s, ' 
            
            %save minima locations.
            
            %write the offsets to the output file.
            if iscell(expt) == 1
                expt = expt{1};
            end
            if strcmpi(to_do, 'profile track') == 1
                
                outfile_name = strcat(expt,'_profiles.txt');

                %make headers
                headers.file_name = file_name;
                headers.run_name  = run_name;
                headers.caller    = [mfilename('fullpath'),'.m'];
                headers.Xmin      = x_im(1);
                headers.Xmax      = x_im(2);

                status = WritePositionChange(outfile_name, 'OpenNew');
                status = WritePositionChange(outfile_name, 'profile header', headers);
                status = WritePositionChange(outfile_name, 'values', prof_mins, {out_ref{:}}, timestamps);
                status = WritePositionChange(outfile_name, 'close');
            end
            
            
            outfile_name = strcat(run_name,'_profiles.mat');
            save(outfile_name, 'profiles_out', 'x_im', 'y_im', 'out_ref', 'timestamps')
            
%             
%             XLS_out = cell(length(uids)+5,length(x_mins)+2);
%             XLS_out(1:3,1:3) =  [{'Input file:'}, {file_name}, {''}; {'Band averaged over:'}, {'X_min'}, num2cell(x_im(1)); {'-'}, {'X_max'}, num2cell(x_im(2))];
%             XLS_out(5,1:2) = [{'File'}, {'Timestamp'}];
%             XLS_out(6:end,:) = [uids num2cell([timestamps-timestamps(1) prof_mins])];
%             %XLS_out(6:end,:) = [uids num2cell([timestamps-timestamps(1) prof_mins+y_mins(1)])];
%            
%             %FIX ME. The xlswrite function seems to change the time stamps
%             %when it writes the files.
%             if iscell(expt) == 1
%                 expt = expt{1};
%             end
%             xlswrite( (strcat(expt,'_profiles.xls')), XLS_out);
%             

        elseif strcmpi(to_do, 'scale') == 1
            close(h)
            ImageScale(prof')
            return
        
        elseif strcmpi(to_do, 'movie') == 1
            close(writerObj);
%             aviobj = close(aviobj);
        end


    case 'disp'

        %%% checks to see if the boxes are selected - if not then makes them
        if exist(box_file_name,'file') == 0
            warning('ImageAnalysis:NoBoxes','No boxes have been selected or the file is not recognised. New boxes are being created.');
            to_run = strcat(mfilename,'(file_name,''boxes'')');
            eval(to_run);
            close all;
        end

        %loads the box positions.
        Boxes = load(box_file_name);
        if isfield(Boxes, 'Image_rotation') == 1
            Rotation = Boxes.Image_rotation;
        else
            Rotation = 0;
        end
        
        %error checking on loaded boxes. Checks to see if the boxes are not
        %integers or overlap -- catching errors generated by old veriosns of this code 
        if round(Boxes.boxX) ~= Boxes.boxX %| round(Boxes.boxY) ~= Boxes.boxY
            error('Horizontal box positions are not integers!')
        elseif size(Boxes.boxX,1) ~= 1 && Boxes.boxX(1,2) == Boxes.boxX(2,1)
            error('Boxes overlap - need updating')
        end
        
        
        % Checks the analysis type and sets where displacement analysis starts.
        % 2 for 'seq' because images referenced to previous image; 1 for 'ref' because all images are analysed relative to middle image
        % This goes here before the matlabpool/parpool is opened so that smpd.exe
        %     does not inherit the opened position_change file - which is then
        %     undeletable until matlab is closed. 
        
        %determines if parallel processing is needed or not.
        if strcmpi(AnalysisVariables.analysis_type,'seq')
%             start = 2;
%             summed_displacements = zeros(1,size(Boxes.boxY,1));
            block_frames = 50; % reads images in blocks to speed to running up.
            parforArg = 1;
            
        elseif strcmpi(AnalysisVariables.analysis_type,'ref') && length(Boxes.boxX) >= 5
%             start = 1;

            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(poolobj)
                poolsize = 0;
            else
                poolsize = poolobj.NumWorkers
            end
            
            % undocumented function to get number of logical core in processor.
            % https://undocumentedmatlab.com/blog/a-few-parfor-tips
            numCores = feature('numcores'); 
            
            parforArg = numCores;
            
            if poolsize == 0
                parpool(parforArg)
            end
            %if matlabpool('size') == 0
            %    matlabpool open
            %end
            
            block_frames = 100;
        else % analysis type = ref and 5 or fewer boxes
%             start = 1;
            block_frames = 100;
            parforArg = 0;
        end
               
        %create array to write in both positionchange headers and SSD files
        headers.file_name = file_name;
        headers.run_name  = run_name;
        headers.caller    = [mfilename('fullpath'),'.m'];
        headers.variables = AnalysisVariables;
        headers.boxes     = Boxes;
        headers.refid     = ref_id;
                
        %opens the displacement output file
        %   -- if it doesn't exist it creates the file and writes the header in it.
        if exist(outfile_name,'file') == 0 && save_disp == 1
            status = WritePositionChange(outfile_name, 'OpenNew');
            
%             disp('no position_change file -- creating header.')
            
            status = WritePositionChange(outfile_name, 'header', headers);
            
        elseif save_disp == 1
            status = WritePositionChange(outfile_name, 'Open');
        end
                
        % For benchmarking
        %         max_frames = 300; tic;
        
        %Dont know why but parfor seems to need this initialisation
        AVmin_type = AnalysisVariables.min_type;
        xBox = round(Boxes.boxX);
        yBox = round(Boxes.boxY); %round the boxes because they might not be integers
%         AVsearch = AnalysisVariables.search;
%         AVspot_removal = AnalysisVariables.spot_removal;
%         AVnan_bg = AnalysisVariables.NaN_bg;
        
        %rotate image if required
        if Rotation ~= 0
            image1 = rotate_images(image1, Rotation);% Im_procs('rotate',image1, Rotation);
        end
        
        %%%pass image to displacement anaylsis for preprocessing -- NaN bright spots and NaN background
        %%%image1 = displacement_analysis('preprocess', image1, [], xBox, yBox, AVsearch, AVspot_removal, AVnan_bg);
        
        %setup memory space for variables.
        %FIXME. Need to set all the used arrays here. 
        offsets = zeros(block_frames, size(yBox,1));
        if strcmpi(AnalysisVariables.analysis_type,'seq')
            summed_displacements = zeros(1,size(Boxes.boxY,1));
        end
        
        disp_args = {};
        disp_args = [disp_args, 'search', AnalysisVariables.search];
        if plots_on == 1
            disp_args = [disp_args, 'plot'];
        end
        if strcmpi(AnalysisVariables.spot_removal,'yes')% == 1
            disp_args = [disp_args, 'spot removal'];
        end
        if strcmpi(AnalysisVariables.NaN_bg,'yes')% == 1
            disp_args = [disp_args, 'nan bg'];
        end
        if strcmpi(AnalysisVariables.offset,'yes')% == 1
            disp_args = [disp_args, 'x offset', 'find'];
        end
        
        %big loop for the displacement analysis
        for j = start : block_frames : max_frames
            
            if j+block_frames > max_frames
                block_frames = max_frames - j + 1;
            end
            
            fprintf(1, '%s; Reading images %4.0f to %4.0f \n', run_name, j, j+block_frames-1);
            second_frames = Im_procs('getimage',file_id, AnalysisVariables.expt_location, j, block_frames, image_size);
            
            % remove dark field if required. Needs to be done before
            % anything else is done because the dark field is direclt from
            % unprocessed images.
            if Subtract_DF == 1
                second_frames = second_frames-repmat(DF_image,[1,1,size(second_frames,3)]);
            elseif Subtract_DF == 2
                second_frames = FFTImageFilter(second_frames, DF_image);
            end
            
            if Rotation ~= 0
                second_frames = rotate_images(second_frames, Rotation);
            end
                
            nos = (1:block_frames) + j - 1;
            
            switch AnalysisVariables.analysis_type
                case 'ref' %use reference image for analysis
                  
                    %pass image to displacement anaylsis for preprocessing -- NaN bright spots and NaN background
                    if j == start
                        image1 = displacement_analysis('preprocess', image1, [], xBox, yBox, disp_args{:});
                    end
                    %this is done here so that the pre-processing for the
                    %displacement experiments can be done in the correct
                    %place as well. 
                    % only needs to be done once. 

                    %Call function to do the work
                    % separate the calls for paralell and linear for loops
                    % to speed up code when small number of boxes
                    if parforArg > 1
                        % This should be replaced with calls to parfeval
                        % because the execution is asynchronous.
                        parfor (x = 1 : block_frames,parforArg) 
                            fprintf(1, '%s; Doing cross-corellation of images %4.0f and %4.0f \n', run_name, ref_id, nos(x));
                            image2 = second_frames(:,:,x);
                            image2 = displacement_analysis('preprocess', image2, [], xBox, yBox, disp_args{:});
                            [offsets(x,:), SSD_temp_array(:,:,x), pos_temp_array(:,:,x)]  = displacement_analysis(AVmin_type, image1, image2, xBox, yBox, disp_args{:}); %<-disp analysis row->% this is needed to find right row for disp analysis version
                        end
                    else
                        for x = 1 : block_frames
                            fprintf(1, '%s; Doing cross-corellation of images %4.0f and %4.0f \n', run_name, ref_id, nos(x));
                            image2 = second_frames(:,:,x);
                            image2 = displacement_analysis('preprocess', image2, [], xBox, yBox, disp_args{:});
                            [offsets(x,:), SSD_temp_array(:,:,x), pos_temp_array(:,:,x)]  = displacement_analysis(AVmin_type, image1, image2, xBox, yBox, disp_args{:}); %<-disp analysis row->% this is needed to find right row for disp analysis version
                        end
                    end  
                    for x = 1 : block_frames %FIX ME: it should be possible to make these arrays without a loop. (e.g. like the time stamp array)
                        out_ref{x} = uids(ref_id); %N.B. out_ref and out_nos need to be cell arrays inorder to cope with any change in length of the strings.
                        out_nos{x} = uids(nos(x));
                    end 
                    out_times = timestamps(nos(1:block_frames));
                    
                case 'seq' %sequential image analysis
                        
					%preprocess the first image
                    image1 = displacement_analysis('preprocess', image1, [], xBox, yBox, disp_args{:});
                            
                    for x = 1 : block_frames
                        
                
                        %get the second image
                        image2 = second_frames(:,:,x);
                        image2 = displacement_analysis('preprocess', image2, [], xBox, yBox, disp_args{:});
                        
                        fprintf(1, '%s; Doing cross-corellation of images %4.0f and %4.0f \n', run_name, ref_id, nos(x));
                        
                        %Call function to do the work
                        [offsets(x,:), SSD_temp_array(:,:,x), pos_temp_array(:,:,x)] = displacement_analysis(AVmin_type, image1, image2, xBox, yBox, disp_args{:}); %<-disp analysis row->% this is needed to find right row for disp analysis version
                        
                        %makes record of files compared for writing to the output file
                        out_ref{x} = uids(ref_id); %N.B. out_ref and out_nos need to be cell arrays inorder to cope with any change in length of the strings.
                        out_nos{x} = uids(nos(x));
                        out_times(x) = timestamps(nos(x));
                        
                        % plot image of moving foils with boxes overlain
                        figure(1)
                        subplot(1,2,1); imagesc(image1); colormap(gray); title(['Image ',num2str(ref_id),' in sequence']);
                        set(gca, 'clim', [0 256])
                        for n = 1 : size(yBox,1)
                            bx_x = xBox(n,1);
                            bx_y = yBox(n,1);
                            bx_h = xBox(n,2)-xBox(n,1);
                            bx_w = yBox(n,2)-yBox(n,1);
                            rectangle('position',[bx_x,bx_y,bx_h,bx_w],'EdgeColor','r');
                        end
                        subplot(1,2,2); imagesc(image2); set(gca, 'clim', [0 256]), title(['Image ',num2str(nos(x)),' in sequence']);
                        for n = 1 : size(yBox,1)
                            bx_x = xBox(n,1);
                            bx_y = yBox(n,1)+offsets(x,n);
                            bx_h = xBox(n,2)-xBox(n,1);
                            bx_w = yBox(n,2)-yBox(n,1);
                            rectangle('position',[bx_x,bx_y,bx_h,bx_w],'EdgeColor','r');
                        end
%                         if nos(x) > 90 && nos(x) < 100
%                             waitforbuttonpress
%                         end

                        %copies the previous image into the 'reference' image position.
                        image1 = image2;
                        ref_id = ref_id+1;

                        %moves the boxes to keep up with changing position of foils in 'first' image
                        summed_displacements = summed_displacements + offsets(x,:);
%                         boxY = boxY0 + [round(summed_displacements); round(summed_displacements)]';
                        yBox = yBox + [offsets(x,:); offsets(x,:)]';


                    end
            end  %end analysis type switch
        
            %write the offsets to the output file. 
            if save_disp == 1
                status = WritePositionChange(outfile_name, 'values', offsets(1:block_frames, :), out_ref(1:block_frames), out_nos(1:block_frames), out_times(1:block_frames));
            end
                     
            %store the SSD values.
            if save_SSD == 1 
                num_lines = block_frames;
                SSD_array(:,:,j-start+1:j-start+num_lines) = SSD_temp_array(:,:,1:num_lines);
                pos_array(:,:,j-start+1:j-start+num_lines) = pos_temp_array(:,:,1:num_lines);
                image_time(j-start+1:j-start+num_lines) = out_times(1:num_lines);
            end
            
        end % block frames loop
        
        %save the SSD values        
        if save_SSD == 1
%             ImageAnalysis_SSDwrite(SSD_file_name, SSD_array, pos_array, image_time, headers)
            ImageAnalysis_SSDio('write', SSD_file_name, SSD_array, pos_array, image_time, headers)
        end
        
        if save_disp == 1
            WritePositionChange(outfile_name, 'close');
        end
        
    otherwise %to_do switch - which process to undertake.
        error('Requested process is not recognised')
end

Im_procs('close',file_id);

fclose all;

toc
end

%% subfunctions

function make_tiffs(image1,name,k)

name = make_im_names(name,k, 'tif');
if 1%max(image1)> 2^8
    image1 = uint8(image1);
else
    image1 = uint16(image1);
end
imwrite(image1, name, 'Compression', 'none')
%asdf = imread(name);
%imagesc(asdf)
%keyboard
end

function make_pngs(image1,name,k)

name = make_im_names(name,k,'png');
if 1%max(image1)> 2^8
    image1 = uint8(image1);
else
    image1 = uint16(image1);
end


imwrite(image1, name, 'Compression', 'none')

end


function out = make_im_names(name,k,ending)
out = strcat(name,'_image',sprintf("%05i",k),'.',ending);
end

function out = rotate_images(images, angle)
%         ims = varargin{1};
%         angle = varargin{2};
        num_ims = size(images,3);
        for x = 1:num_ims
            out(:,:,x) = imrotate(images(:,:,x),angle,'bilinear');
        end
%         varargout{1} = rot_ims;
end %rotate_images



function out=polymin(data)

range = 10;
% data

% data(data>mean(data))=[];
loc = find(min(data) == data);
loc = loc-range:loc+range;

loc(find(loc <= 0)) = [];
loc(find(loc > length(data))) = [];

dat = data(loc);

p  = polyfit(1:length(dat),dat,3);
diff_p = [p(1)*3 p(2)*2 p(3)];
p_roots = roots(diff_p);

p_roots = real(p_roots);

p_roots(p_roots>length(dat)) = [];
p_roots(p_roots<1) = [];


if length(p_roots)>=1
   valus = polyval(p,p_roots);
   
   [~,b] = min(valus);
   out = p_roots(b);
else
    arr = 0:0.01:length(dat);
    p_calc = polyval(p,arr);
    out = arr(find(min(p_calc)==p_calc))+loc(1)-1;
    %out = p_roots;
end

out = out(1)+loc(1);

if isempty(out)
    plot(data)
    disp(out)
    keyboard
end
%  figure(3), clf
% plot(data), hold on,
% plot(out,min(data),'rx')
% pause
% keyboard




end

%% VERSION info:
%  - v2.4 -- 21st September 2019
%       - Added fft filter to the images. As a background removal option.
%       Calls FFTImageFilter.m
%  - v2.3.1 -- 23rd April 2019
%       - fixed error in passing arguments to displacement_analysis. 
%  - v2.3 -- 10th April 2019
%		- unpdated script to use displacement_analysis version 4. Which as a different syntax to version 3 and earlier.
%		- also made (minor) structural changes to make code more efficient.
%  - v2.2 -- 10th April 2019
%       - replaced matlabpool with parpool commands. In so split 'ref' call
%       in 'disp' loop into parallel and linear options.
%  - v2.1 -- 19th August 2017
%       - Added option to subtract dark field from images before processing.
%  - v2.0.1 -- 19th May 2017
%       - bug fixes for part of script in linux
%  - v2.0.1 -- 11th January 2017
%       - edit parts of the movie options so that it works. Replaced avifile with VideoWriter.
%  - v2.0 -- 20th June 2016
%       - Move image specific functions into their own separate files. 
%       - Separated location specific image manipulations into their own files.
%       - Move position_change writing code into new script called WritePositionChange.m.
%       - Added lines so that the default varaibles can be changed -- for runnning loops over all reference images. 
%       - Separated out SSD file save into a separate file.
%       - Added switches at the top for saving of output files, rather than having to do a string compare each time.
%  - v1.9.1 -- 9th May 2016
%       - Made all the variables feeding into the displcement_analysis.m consistent. 
%  - v1.9 -- 26 April 2016 onwards...
%       - Tidy up and complete documentation. As part of this deleted all the unneeded developlemt versions of the script. 
%       - Wrote file documenting how the sctipt works.
%  - v1.8.2 -- 27th Jan 2016
%       - Chnaged error checking on input boxes to allow updated
%       SelectBoxes which now returns non integer box positions. 
%       - Round boxes.
%  - v1.8.1 -- 30th November 2015
%       - Fixed error which made the SSD output array replicate data if the number of images was not a multiple of block_frames.
%       - Changed writing of 'position_change' file so that the directory name is not writted to the output file with the image name ,
%       this does not change the netcdf output only the tiff and other single image output.  
%  - v1.8 -- 8th October 2015
%       - Added the ability to rotate the images (if required) and process them.
%  - v1.7.1 -- 17 Sept 2014
%       - Tidied up code a bit and made load load variables to a structure rather than individual variables.
%  - v1.7 -- 2 Sept 2014
%       - Added ability to save SSD calculated along with displacements. The options for saving it are all 
%       in this script and outputs from Dislplacement_analysis are fixed.
%  - v1.6.7 -- 12 August 2014
%       - changes to make script work on MacOs and linux. Generally these changes are
%       to how I address file names. The script now uses filesep rather than a hard coded '/'
%       and fullfilename as a function to make unique file addresses.
%  - v1.6.6 -- 11 July 2014
%       - changed parallel options so that only parallel if there are a
%       large number of boxes. Also adjusted block frames to correspond to
%       the efficiencies of the parrarel or not analysis.
%  - v1.6.5 -- 13 Jan 2014
%       - added option for finding scale of images -- needs a series of images and uses ImageScale.m
%  - v1.6.4 -- 9 Jan 2014
%       - minor changes so that can accept fullstops in file names that are not the pre file extension fullstop.
%  - v1.6.3 -- 27 Nov 2013
%       - Changed profile option so that it plots an profile of the average intensity across a selected area. 
%  - v1.6.2 -- 21 Feb 2013
%       - Added couple of lines to make the maximum value in the images 2^8. This is mostly for the sidestation camera at the NSLS but
%       should be made into a separate option.
%  - v1.6.1 -- 21 Feb 2013
%       - Added code to rotate APS images into correct orientation
%  - v1.6 -- 13 Jan 2013
%       - Added code to parse files because it can take a long time to read the files when they are first opened.
%  - v1.5b -- 8th January 2013
%       - Added code to make profile through band of image sequence	
%       - Added code to deal with files without the NSLS naming system
%       - Added code to pre-process DESY images
%       - Changed movie code so that it gets the frames per second from the time stamps and removed .nc. from movie file name
%  - v1.5a -- 16th Nov 2012
%       - Added basic support for 'movie' run producing an avi file of the
%       data. Note the nasty hack to avoid having to use experiment files
%       from conductivity experiments. Also, we should get the frame rate
%       from somewhere so this plays in real time and remove the .nc. from
%       the middle of the filename.
%  - v1.5 -- 1st July 2012
%       - Added code for parallel processing of images when using a reference image. As part of this separated the sequential and
%       reference analysis options so now have figure for sequential analysis.
%       - Forced out_nos and out_ref to be cell arrays. This way they can cope with changes in length of the strings fed to them.
%  - v1.4h -- Feb 2012
%       - added check that the BoxX and boxY numbers are integers and the boxes do not overlap.
%  - v1.4g -- 25 January 2012
%       - added preprocess loop for displacement analysis. Sorts image1 for processing so not repeated each time.
%            Needs displacement_analysis version 3.8 or later.
%  - v1.4f -- 5 December 2011
%       - changed function name for selecting boxes to SelectBoxes
%  - v1.4e -- 24 November 2011
%       - replaced a '\r' with '\n' in the header creation routine. This was making problems in Phases.m
%  - v1.4d -- 12 October 2011
%       - Added line forcing images to be double precision numbers
%  - v1.4c -- 13 Feb 2011
%       - Added options to feed spot removal to displacement_analysis.m
%  - v1.4b -- 30 Jan 2011
%       - rewrote netcdf_func/manipualte image using 'typecast' which makes it much faster
%  - v1.4a -- 29 Jan 2011
%       - fixed error which had script looking for exptfile in wrong place.
%  - v1.4
%       - Moved header creation for output file from displacement_analysis (v3.5) into this file.
%       - Added new lines to header which record the details of the analysis options
%       - Moved how far search into setting of ImageAnalysis
%       - Analysis variables are now stort in 'expt-name'.mat and read from there.
%       - Removed all superflous lines of code.
%  - v1.3 - 1st numbered version 1 Jan 2011
%     - Added options for sequential and referential analysis (October 2010).
%     - Change so that script can take different file formats easily (August 2010).
%           For a sequence of tiff images this requires a list file as well.
%
%  - prior to July 2010.
%     - made to work with matlab native NetCDF routines.