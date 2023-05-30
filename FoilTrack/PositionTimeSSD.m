% PositonTimeSSD
%   Calculate the position and displacement rate of the foils using SSD data produced by ImageAnalysis script. 
%   It is designed to be used on data which is too noisy to do on a per image
%   basis or for data which has very small smplitudes. 
%
%   syntax: PositonTimeSSD(SSD_file_name.mat, options)
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

function PositionTimeSSD(file_name, varargin)

tic
fprintf(1, '\n Dispalcements in %s \n', file_name);

warning off all

%% get settings and variables.

%set default options
window      = 6;
m           = 6;
n           = 3;
begin       = -1;
cease       = -1;
plots_on    = 0;
output_save = 1;
overwrite  = 0;

% prase varargin...
iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'window'
            window = varargin{iarg+1};
            iarg = iarg + 2;
        case 'surface coefficients'
            m = varargin{iarg+1};
            n = varargin{iarg+2};
            iarg = iarg + 3;
        case 'begin'
            begin = varargin{iarg+1};
            iarg = iarg + 2;
        case 'cease'
            cease = varargin{iarg+1};
            iarg = iarg + 2;
        case 'search'
            search = varargin{iarg+1};
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
[~, root, ~] = fileparts(file_name);
root_length = strfind(root,'_SSD');
root = root(1:root_length-1);
tail = '_SSD_displacements_rates';
out_file_name1 = [root,'_SSDwindow_positions.txt']; %removes 'position_change' from the end of the input data file name.
out_file_name2 = [root,'_SSDwindow_lengths.txt']; %removes 'position_change' from the end of the input data file name.
out_file_name3 = [root,'_SSDwindow_PositionChangeRate.txt']; %removes 'position_change' from the end of the input data file name.
out_file_name4 = [root,'_SSDwindow_LengthChangeRate.txt']; %removes 'position_change' from the end of the input data file name.


%% data setup 

% % define order of polynomial surface
if ~isfinite(n) %check to see if n and m were defined with the inputs.
    m = 6; %polynomial order in displacement
    n = 2; %order in time direction %FIX ME THIS HAS BEEN CHANGED FROM 3.
end

% read data file.
% data = load(file_name, 'pos_array', 'SSD_array', 'image_time', 'boxX', 'boxY');
data = ImageAnalysis_SSDio('read', file_name, 'pos_array', 'SSD_array', 'image_time', 'boxX', 'boxY');

% extract data and reorient into 'correct' orientation
pos_array = permute(data.pos_array, [2 3 1]); %pos_array is absolute position values in image
SSD_array = permute(data.SSD_array, [2 3 1]);
start_time = data.image_time(1);
image_time = data.image_time(:) - start_time; %make time stamps start from zero
%get size properties of the data arrays.
dat_size = size(pos_array);
if length(dat_size) == 2
    number_boxes = 1;
else
    number_boxes = dat_size(3);
end
% 
% %get box reference positions
% if strcmpi(file_name(end-2:end),"SSD")
%     boxes = load([file_name(1:end-4), '_boxes']);
% else
%     boxes = load([file_name, '_boxes']);
% end
% data.boxX = boxes.boxX;
% data.boxY = boxes.boxY;



%move SSD loacations be around the reference position. (i.e. -n to n)
offset_values = pos_array((dat_size(1)+1)/2, round(dat_size(2)/2), :);
pos_array = pos_array - repmat(offset_values,[dat_size(1), dat_size(2), 1]);

% Perform data pre-processing unique to experiment.
% the file 'Experiment_data_processing.m' needs to be in the current directory.
% It contains processing unique to the experiment e.g. for the Zn_? series of experiments the cutting of the begining and end of the data set is done here. 
if exist('Experiment_data_processing.m','file')
    preprocfunc = [mfilename,'_DataCut'];
    [pos_array, SSD_array, image_time] = Experiment_data_processing(preprocfunc,pos_array, SSD_array, image_time);

    %If too much has been cut exit otherwise recalcualte the dimension arrays, just incase they have changed.
    if numel(pos_array) == 0
        fprintf('No data left after data size reduction. Exiting \n\n')
        return
    else
        dat_size = size(pos_array);
        number_boxes = dat_size(3);
    end
else
    fprintf('Experiment_data_processing does not exist on path. No modifications to the data set.\n')
end

% Smooth the time stamps (if necessary) and make times start from 0
arg_to_pass = {}; if plots_on == 1, arg_to_pass = {arg_to_pass{:}, 'plot'}; end
image_time = TimeStampSmooth(image_time, arg_to_pass{:}); %required for tiff image based data sets. If netcdf based the returns same array as input.
image_time = repmat(image_time',dat_size(1),1); %makes image_time array equivalent size to other arrays


%debugging function.
if (0)
    fprintf('Size of data array      : %d %d %d \n', size(pos_array,1), size(pos_array,2), size(pos_array,3))
    fprintf('Size of time stamp array: %d %d    \n', size(image_time,1), size(image_time,2))
    keyboard;return
end



%% SSD shape for individual boxes
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
    
    
    pos_now  = pos_array(:,j:j+window-1,:);
    SSD_now  = SSD_array(:,j:j+window-1,:);
    time_now = image_time(:,j:j+window-1);


    %fit the polynomial shape without accounting for sinusoidal vairations with
    %time.
    
    % generate polynomial coefficients for foils.
    A = zeros(1,polyvalnm_ncoef(m,n));
    
    t_mean(j) = mean(time_now(:));
    X = time_now-t_mean(j);
    X = X(:);
    
    opts = statset('TolX', 1e-12);
    
    for x = 1 : number_boxes
        
        %print running tally of where we are. 
        str = sprintf('  Box %i of %i pair %i', x, number_boxes, j);
        del = sprintf(repmat('\b',1,length(str)));
        fprintf(str)
        
        Y = pos_now(:,:,x);
        Y = Y(:);
        in = [X Y];
        
        Z = SSD_now(:,:,x);
        Z = Z(:);
        
        model = NonLinearModel.fit(in, Z, @poly_model, A, 'Options', opts);
        A_all(x,:) = model2array(model.Coefficients(:,1));
    
        % Assuming that the foils are all about the same shape use the coefficients for the previous foil to fit the next one.
        A = A_all(x,:); 
    
        %if required plot data
        if 0%plots_on == 1
            surface_plot(time_now, pos_now(:,:,x), SSD_now(:,:,x), model.Fitted, NaN)
            pause
        end
        fits(j,x,:) = model.Fitted;
        surf_coef_all(j,x,:) = model.Coefficients.Estimate;
        
        fprintf(del)
    
        %calculate position and gradient
        
        %differentiate the surface fit
        surf_coef = polyvalnm_coef2mat(surf_coef_all(j,x,:), n,m);
        surf_diff = polyvalnm_diff(surf_coef, 2,1);
        
        %solve for position at average time
        mean_time_now = mean(X);

        displacement(j,x,:) = polyvalnm_solve2(surf_diff, 0, mean_time_now);


        %% calculate the displacement rate.
        multi_displacements = polyvalnm_solve2(surf_diff, 0, min(X):(max(X)-min(X))/200:max(X));
        min_poly = polyfit(min(X):(max(X)-min(X))/200:max(X),multi_displacements,n);

        slope = polyder(min_poly);
        slope_at_mean(j,x,:) = polyval(slope,mean_time_now);

        if 1%plots_on==1

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
    timestamp(j) = t_mean(j);
        
    
    %if required plot data
    if 0%plots_on == 1
        surface_plot2(X, pos_now, SSD_now, squeeze(fits(j,:,:)), NaN)
        pause
    end
end

figure
plot(slope_at_mean)


%% combine data and prepare for output
position = repmat(mean(data.boxY,2)',size(displacement,1),1) + displacement;
l  = position(:,2:end)-position(:,1:end-1);

length_change_rate = slope_at_mean(:,2:end)-slope_at_mean(:,1:end-1);

figure
subplot(1,2,1)
plot(l)
subplot(1,2,2)
plot(length_change_rate)



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
AnalysisVariables.surface_coef = [m,n];

%get the box positions
box_file_name = [root,'_boxes.mat']; %the file in which the box positions are stored
%loads the box positions.
Boxes = load(box_file_name);
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
headers.refid     = data.refid;
        

%% write output files
for x = 1:4

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
        v = slope_at_mean;

    elseif x==4
        % length change rate
        of = out_file_name4;
        v = length_change_rate;

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
        status = WritePositionChange(of, 'values', v(j,:), refs(j), image_ids(j), timestamp(j)+start_time);
    end
    WritePositionChange(of, 'close');

end





% % %opens the position output file
% % %   -- if it doesn't exist it creates the file and writes the header in it.
% % if exist(out_file_name1,'file') == 0
% %     status = WritePositionChange(out_file_name1, 'OpenNew');
% % else
% %     status = WritePositionChange(out_file_name1, 'Open');
% % end
% % status = WritePositionChange(out_file_name1, 'header', headers);
% % 
% % %write the positions
% % for j = process_start : processs_end-window+1
% %     %write the offsets to the output file. 
% %     status = WritePositionChange(out_file_name1, 'values', position(j,:), refs(j), image_ids(j), timestamp(j));       
% % end
% % WritePositionChange(out_file_name1, 'close');
% % 
% % 
% % %opens the length output file
% % %   -- if it doesn't exist it creates the file and writes the header in it.
% % if exist(out_file_name2,'file') == 0
% %     status = WritePositionChange(out_file_name2, 'OpenNew');
% % else
% %     status = WritePositionChange(out_file_name2, 'Open');
% % end
% % headers.boxes.boxY = Boxes.boxY(2:end,:)-Boxes.boxY(1:end-1,:);
% % headers.boxes.boxX = Boxes.boxX(1:end-1,:);
% % status = WritePositionChange(out_file_name2, 'header', headers);
% % 
% % %write the positions
% % for j = process_start : processs_end-window+1
% %     %write the offsets to the output file. 
% %     status = WritePositionChange(out_file_name2, 'values', l(j,:), refs(j), image_ids(j), timestamp(j));       
% % end
% % WritePositionChange(out_file_name2, 'close');
% % 
% % 
% % 
% % %write the length change rate output file
% % %   -- if it doesn't exist it creates the file and writes the header in it.
% % if exist(out_file_name3,'file') == 0
% %     status = WritePositionChange(out_file_name3, 'OpenNew');
% % else
% %     status = WritePositionChange(out_file_name3, 'Open');
% % end
% % headers.boxes.boxY = Boxes.boxY(2:end,:)-Boxes.boxY(1:end-1,:);
% % headers.boxes.boxX = Boxes.boxX(1:end-1,:);
% % status = WritePositionChange(out_file_name3, 'header', headers);
% % 
% % %write the positions
% % for j = process_start : processs_end-window+1
% %     %write the offsets to the output file. 
% %     status = WritePositionChange(out_file_name3, 'values', l(j,:), refs(j), image_ids(j), timestamp(j));       
% % end
% % WritePositionChange(out_file_name3, 'close');


%%%%%% end of function. Subfunctions here after.

%% subfunctions


    function [out] = poly_model(coef, predictors)
        
        t = predictors(:,1);
        displ = predictors(:,2);
        
        %reshape polynomial coefficients
        coef = polyvalnm_coef2mat(coef,n,m);
        
        %caclulate polynomial surface using adjusted x values.
        out = polyvalnm(coef,displ,t); 
               
    end %poly_model


end %Phases_SSD


function out = model2array(in)
% this function covers for differences in the output of the fitting
% routines between matlab 2012a and later versions.

if isa(in, 'table') == 1
    out = table2array(in);
% elseif isa(in, 'dataset') == 1 
%     out = cell2mat(dataset2cell(in));
else
    out = double(in);
end

end

function surface_plot(time, posit, SSDs, fit, per )

view_ogn = [-68,25];

fg = figure(2);

%set(fg,'Renderer','OpenGL');

%plot data
ax(1) = subplot(1,3,1);
    %surf(time, posit,SSDs), hold on
    %scatter3(time(:), posit(:),SSDs(:), 1, 'k'), hold off
    scatter3(time(:), posit(:),SSDs(:), 4, SSDs(:)), hold off
    view( view_ogn );
    shading interp
    ylabel({'Offset (pixels)'}, 'Rotation', -10)
    xlabel('Time (s)', 'Rotation', 60)
    cb(1) = colorbar('horizontal');
    set(get(cb(1),'ylabel'),'String','Intensity^2');
    zlabel('SSD (observations, relative to reference image)')
    title('a. SSD data')

    
    set(ax(1), 'XTickLabelRotation', -10, 'YTickLabelRotation', 0, 'ZTickLabelRotation', -10, 'FontSize', 12)
    
%     try
%         h = rotate3d;
%         set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
%         set(gcf, 'ResizeFcn', @align_axislabel)
%         align_axislabel([], gca)
%         %axislabel_translation_slider;
%     end
%plot model
fit = reshape(fit,size(time));

ax(2) = subplot(1,3,2);
    surf(time, posit, fit)
    view( view_ogn );
    shading interp
    ylabel({'Offset (pixels)'}, 'Rotation', -10)
    xlabel('Time (s)', 'Rotation', 60)
    zlabel('Model SSD')
    title('b. Model fit to data')
    cb(2) = colorbar('horizontal');
    set(get(cb(2),'ylabel'),'String','Intensity^2');
    set(ax(2), 'XTickLabelRotation', -10, 'YTickLabelRotation', 0, 'ZTickLabelRotation', -10, 'FontSize', 12)
%     try
%         h = rotate3d;
%         set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
%         set(gcf, 'ResizeFcn', @align_axislabel)
%         align_axislabel([], gca)
%         %axislabel_translation_slider;
%     end 
    
%plot Residuals      
Resids = (SSDs - fit);%./(SSDs);
%Resids = (sqrt(SSDs) - sqrt(fit))./sqrt(SSDs);

ax(3) = subplot(1,3,3);
	d = time(:);
    e = posit; e = e(:);
    f = Resids(:);
    if 1;%isfinite(per) == 0
        %plot residuals as function of time.
        scatter3(d, e, f, 12, f) %def because has to be plotted as vectors
        xlabel('Time (s)', 'Rotation', 60)
    else
        %plot residuals as function of phase angle.
        %scatter3(rem(d,per)./(per)*(2*pi), e, f, 12, f) %def because has to be plotted as vectors
        scatter3(d, e, f, 12, f) %def because has to be plotted as vectors
        xlabel({'  Phase  ';'(radians)'}, 'Rotation', 60)
        %xlim([0,2*pi])
    end
    view( view_ogn );
    ylabel({'Offset (pixels)'}, 'Rotation', -10)
    zlabel('Residuals (SSD - model)')
    title('c. Residuals')
    colormap(ax(3),'cool')
    cb(3) = colorbar('horizontal');
    set(get(cb(3),'ylabel'),'String','Intensity^2');
    set(ax(3), 'XTickLabelRotation', -10, 'YTickLabelRotation', 0, 'ZTickLabelRotation', -10, 'FontSize', 12)
%     try
%         h = rotate3d;
%         set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
%         set(gcf, 'ResizeFcn', @align_axislabel)
%         align_axislabel([], gca)
%         %axislabel_translation_slider;
%     end
    
%link all the subplots together.
%FIX ME. for some reason this does not work. Something to so with using both surf and scatter3.
   %      Link = linkprop(ax, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
   %      setappdata(fg, 'StoreTheLink', Link);

%keyboard%pause(.2)

end


function surface_plot2(time, posit, SSDs, fits, per )

view_ogn = [-68,25];

number_foils = size(SSDs,3);

fg = figure(2);

mainfig = figure(2);
tabgroup = uitabgroup(mainfig, 'Position', [.01 .01 .98 .98]);

for y = 1:number_foils
    
    tab(y)=uitab(tabgroup,'Title', sprintf('Foil_%i', y));
    axes('parent',tab(y))
                   
    ax(y,1) = subplot(1,3,1);
    a = time;
    b = posit(:,:,y);
    c = SSDs(:,:,y);
    scatter3(a(:), b(:),c(:), 4, c(:)), hold off
    view( view_ogn );
    shading interp
    ylabel({'Offset (pixels)'}, 'Rotation', -10)
    xlabel('Time (s)', 'Rotation', 60)
    cb(1) = colorbar('horizontal');
    set(get(cb(1),'ylabel'),'String','Intensity^2');
    zlabel('SSD (observations, relative to reference image)')
    title('a. SSD data')

    set(ax(y,1), 'XTickLabelRotation', -10, 'YTickLabelRotation', 0, 'ZTickLabelRotation', -10, 'FontSize', 12)
    
%     try
%         h = rotate3d;
%         set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
%         set(gcf, 'ResizeFcn', @align_axislabel)
%         align_axislabel([], gca)
%         %axislabel_translation_slider;
%     end

    %plot model
    if size(time,2)==1
        time = reshape(time,size(posit(:,:,y)));
    end
    fit = reshape(fits(y,:),size(time));

    ax(y,2) = subplot(1,3,2);
    surf(time, posit(:,:,y), fit)
    view( view_ogn );
    shading interp
    ylabel({'Offset (pixels)'}, 'Rotation', -10)
    xlabel('Time (s)', 'Rotation', 60)
    zlabel('Model SSD')
    title('b. Model fit to data')
    cb(2) = colorbar('horizontal');
    set(get(cb(2),'ylabel'),'String','Intensity^2');
    set(ax(y,2), 'XTickLabelRotation', -10, 'YTickLabelRotation', 0, 'ZTickLabelRotation', -10, 'FontSize', 12)
%     try
%         h = rotate3d;
%         set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
%         set(gcf, 'ResizeFcn', @align_axislabel)
%         align_axislabel([], gca)
%         %axislabel_translation_slider;
%     end 
    
    %plot Residuals
    Resids = (SSDs - fit);%./(SSDs);
    %Resids = (sqrt(SSDs) - sqrt(fit))./sqrt(SSDs);
    
    ax(y,3) = subplot(1,3,3);
    d = time(:);
    e = posit(:,:,y); e = e(:);
    f = Resids(:,:,y); f = f(:);
    if 1;%isfinite(per) == 0
        %plot residuals as function of time.
        scatter3(d, e, f, 12, f) %def because has to be plotted as vectors
        xlabel('Time (s)', 'Rotation', 60)
    else
        %plot residuals as function of phase angle.
        %scatter3(rem(d,per)./(per)*(2*pi), e, f, 12, f) %def because has to be plotted as vectors
        scatter3(d, e, f, 12, f) %def because has to be plotted as vectors
        xlabel({'  Phase  ';'(radians)'}, 'Rotation', 60)
        %xlim([0,2*pi])
    end
    view( view_ogn );
    ylabel({'Offset (pixels)'}, 'Rotation', -10)
    zlabel('Residuals (SSD - model)')
    title('c. Residuals')
    colormap(ax(3),'cool')
    cb(3) = colorbar('horizontal');
    set(get(cb(3),'ylabel'),'String','Intensity^2');
    set(ax(y,3), 'XTickLabelRotation', -10, 'YTickLabelRotation', 0, 'ZTickLabelRotation', -10, 'FontSize', 12)
%     try
%         h = rotate3d;
%         set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
%         set(gcf, 'ResizeFcn', @align_axislabel)
%         align_axislabel([], gca)
%         %axislabel_translation_slider;
%     end
    
%link all the subplots together.
%FIX ME. for some reason this does not work. Something to so with using both surf and scatter3.
   %      Link = linkprop(ax, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
   %      setappdata(fg, 'StoreTheLink', Link);

%keyboard%pause(.2)

end
end

function RedChi2 = ReducedChiSquared(obs,model, exclusions)
%copied from chi2test on mathworks website.

obs = obs(~exclusions);
model = model(~exclusions);

%Normalize x
N = length(obs);
x = obs - model;

% x = (x-mean(x))/std(x); %standardization
x = x/std(x);

xp = [-inf, -1.6:.4:1.6, inf]; %tested bins
E = 0.5*N*diff(erfc(-xp/sqrt(2))); %expected frequency
S = histc(x, xp);

O = S(1:end-1); %%observed frequency
% figure(1), plot(xp(2:end),E,'k-',xp(2:end),O,'k.');  pause

Chi2=sum((E-O').^2./E); %statistics
d = (length(xp)-1) - 2; %%degrees of freedom
RedChi2 = Chi2/d; %%reduced Chi^2 value.

end %REdChiSQ

