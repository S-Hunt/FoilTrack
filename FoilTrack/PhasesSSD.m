% PhasesSSD
%   Calculate the phases from SSD data produced by ImageAnalysis script. 
%   It is designed to be used on data which is too noisy to do on a per image
%   basis and for data which has very small smplitudes. 
%
%   syntax: PhasesSSD(SSD_file_name.mat, options)
%     options: 'period', value  -- estimated value for period
%              'fixed'          -- force the period be that given by period.
%              'plot'           -- plot the validation figures.
%              'surface coefficients' -- set the order of the surface
%                                                       coefficients for the fitting 
%
%   See Also: ImageAnalysis, PhaseCheck, Phases

%   $ version: 2.5.1 $ 9th April 2019 $
%       Simon Hunt. 2015 - 2019

% This version is to make sure the correct length of samples is reported in
% the output file. 

function PhasesSSD(file_name, varargin)

tic
fprintf(1, '\n Phases of %s \n', file_name);

warning off all

%% get settings and variables.

%set default options
fixed = 0;
refine = 1;
input_period = NaN;
m = NaN;
n = NaN;
plots_on = 0;
output_save = 1;

% prase varargin...
if nargin == 2 && ~ischar(varargin{1})%isfinite(varargin{1}) == 1 
    %this is a catch for the historical syntax. which expected just a
    %number for the period guess...
    warning on
    warning('Inputting only the a number for the period is old syntax')
    warning off
    input_period = varargin{:};
else
    %parse new syntax.
    iarg = 1;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'period'
                input_period = varargin{iarg+1};
                iarg = iarg + 2;
            case 'surface coefficients'
                m = varargin{iarg+1};
                n = varargin{iarg+2};
                iarg = iarg + 3;
            case 'fixed'
                fixed = 1;
                iarg = iarg + 1;
            case 'no refine'
                refine = 0;
                iarg = iarg + 1;
            case 'no save'
                output_save = 0;
                iarg = iarg + 1;
            case 'plot'
                plots_on = 1;
                iarg = iarg + 1;
                warning on
                warning('The plotting bits have not been tested and are not complete')
                warning off
            otherwise
                error(['Unknown option: ' varargin{iarg}]);
        end
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

% if the variables are not read from 'expt_name'.mat run AnalysisOptions to generate them.
if sum(isfield(vars,filevariables)) ~= (length(filevariables))
    [out] = AnalysisOptions;
    load(varfilename, filevariables{:});
end

%make filename for output file
[~, root, ~] = fileparts(file_name);
root_length = strfind(root,'_SSD');
root = root(1:root_length-1);
tail = '_SSD_sine_fits';
out_file_name = [root,tail]; %removes 'position_change' from the end of the input data file name.
Fits_file_name = [root,'_SSD_FITS.mat'];

%parse period given options
if strcmpi(input_period, 'title')==1
    input_period = expt{4};
    
elseif ~strcmpi(input_period, 'title') && ~isnumeric(input_period)
    error(['Error parsing inputs.\n The period is ''', input_period, ''' is not recognised.']);
    
elseif fixed == 1 && ~isfinite(input_period)
   error('Error parsing inputs.\n The period is fixed but there is no value to fix it to.') 
       
end

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
image_time0 = data.image_time(1);
image_time = data.image_time(:) - image_time0; %make time stamps start from zero

%get size properties of the data arrays.
dat_size = size(pos_array);
if length(dat_size) == 2
    number_boxes = 1;
else
    number_boxes = dat_size(3);
end

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


%% get frequency of data by fft. 
if isfinite(input_period)==0 && fixed == 0
    
    % Create initial guess for period by fft if one is not provided from varargin
    disp('FFT of data for rough period')
    
    %input array for fft
    if plots_on == 1
        arr = {'plotfigs'};
    else
        arr = {};
    end
    
    %Call Phases_fft to calcualte the most likely period.
    [ input_period bigger smaller input_phase ] = Phases_fft( SSD_array, image_time, 'poly', arr{:});
    input_period = input_period *.98;
    %amplitudes are not exctracted because the amplitudes in the SSD do not
    %correspond simply to the displacement of the foils
    
else
    %if not getting the range and phases from fft need to define them.    
    
    %define range of periods.
    bigger = input_period*1.02;
    smaller = input_period*0.98;
    
    %define input phases.
    input_phase = ones(1, number_boxes);
end


%% SSD shape for individual boxes
 disp('Fit polymonial shape for each foil')
%fit the polynomial shape without accounting for sinusoidal vairations with
%time.

% generate polynomial coefficients for foils.
A = zeros(1,polyvalnm_ncoef(m,n));

X = image_time;
X = X(:);

opts = statset('TolX', 1e-12);

for x = 1 : number_boxes
    
    %print running tally of where we are. 
    str = sprintf('  Box %i of %i', x, number_boxes);
    del = sprintf(repmat('\b',1,length(str)));
    fprintf(str)
    
    Y = pos_array(:,:,x);
    Y = Y(:);
    in = [X Y];
    
    Z = SSD_array(:,:,x);
    Z = Z(:);
    
    model = NonLinearModel.fit(in, Z, @poly_model, A, 'Options', opts);
    A_all(x,:) = model2array(model.Coefficients(:,1));

    % Assuming that the foils are all about the same shape use the coefficients for the previous foil to fit the next one.
    A = A_all(x,:); 

%     %if required plot data
%     if plots_on == 1
%         surface_plot(image_time, pos_array(:,:,x), SSD_array(:,:,x), model.Fitted, NaN)
%         pause
%     end
    fits(x,:,:) = model.Fitted;
    
    fprintf(del)
end

%if required plot data
if plots_on == 1
    surface_plot2(image_time, pos_array, SSD_array, fits, NaN)
    pause
end

%% Refine period.
%if the preiod is to be solved for refine the output of the fft or input
%period to get the period for the global fit close to the 'real' value. 

if refine == 1 && fixed == 0 
    disp('Refine period')
    steps = 41;
    stp = (bigger - smaller)/(steps-1);
    Ts = smaller:stp:bigger ;
    
    X = image_time;
    X = X(:);
        
    %pre-make the output arrays.
    refined_phase = zeros(steps,number_boxes);
    refined_amp   = zeros(steps,number_boxes);
    SSD           = zeros(steps,number_boxes);
    
    for y = 1 : steps %loop over the predefined periods.
        
        T = Ts(y);
         
        for x = 1 : number_boxes %loop over the boxes one at a time.
            
            %create input arrays
            Y = pos_array(:,:,x);
            Y = Y(:);
            A = A_all(x,:);
            num_coef = length(A);
            
            in = [X Y];
            in = [in zeros(size(Y))];
            in(1,3) = T;
            in(2,3) = num_coef;
            in(3:3+num_coef-1,3) = A;
            
            Z = SSD_array(:,:,x);
            Z = [Z(:)];
                        
            %fit model to data
            model = NonLinearModel.fit(in, Z, @SSD_fit_fixedT, [input_phase(x), .1], 'Options', opts);
            
            %collate output data
            refined_phase(y,x)       = model2array(model.Coefficients(1,1));
            refined_amp(y,x)         = model2array(model.Coefficients(2,1));
            %SSD(y,x)                 = model.SSE;
            SSD(y,x)                 = model.SSE ./ sum(Z.^2);
            
            %if required plot data -- for debugging 
            if (0)
                surface_plot(image_time, pos_array(:,:,x), SSD_array(:,:,x), model.Fitted, T)
            end
            
        end
    end

    % find the best period by assuming it has the smallest SSD.
    %FIX ME. IF there are lots of foils it should probably be done only using
    %the smallest SSDs
    
    [~, locs] = min(SSD);
    loc_best = round(mean(locs));
    
    % sumSSD = sum(SSD,2);
    % loc_best = find(min(sumSSD) == sumSSD);
%     [~, loc_best]=min(SSD,[],1);
%     loc_best = mode(loc_best);
    
    if plots_on == 1
        figure
        plot(Ts, SSD)
        xlabel('Assumed period (s)')
        ylabel('SSD')
        title('SSD calculated at fixed periods')
        xlim([min(Ts), max(Ts)])
        %keyboard
        pause
    end
    
%     input_period          = Ts(loc_best);
    input_period          = mean(Ts(locs));
    input_phase           = refined_phase(loc_best,:);
    input_amplitude       = refined_amp(loc_best,:);
    SSD_best             = SSD(loc_best,:);
    
else
    
    input_amplitude       = 0.1 * ones(size(input_phase));
end



%% global fit.
% or if the period is fixed fit the individual boxes. 
disp('Global fit of the good data')

%find the cleanset foils on which to perform global fit.
%defines a parameter called to_use which lists the foils to use for this.
if fixed == 1
    to_use = [];
    
elseif number_boxes > 10
    %find foils with smallest SSD
    use = min([number_boxes, 10]);
    [~,order] = sort(SSD_best);
    to_use = sort(order(1:use));
    
else
    to_use = 1 : number_boxes;
    
end
number_used = length(to_use);

%global fit to the boxes to be used.
%If the period is fixed this section is skipped over.
if ~isempty(to_use)
    
    %create input arrays
    %SSD data to be matched
    model_data = SSD_array(:,:,to_use);
    model_data = model_data(:);
    
    % time and displacement arrays. (predictors)
    allY = pos_array(:,:,to_use);
    allX = repmat(image_time, [1 1 length(to_use)]);
    whichbox = permute(repmat(1:length(to_use), [dat_size(1) 1 dat_size(2)]), [1 3 2]);
    in = [allX(:) allY(:) whichbox(:) zeros(size(allX(:)))];
    
    %add other variables to 'in'.
    in(1,4) = length(to_use);
    in(2,4) = m;
    in(3,4) = n;
    
    %find same image comparison data to exlude from the fitting.
    loc_same = (model_data == 0);
    
    %coefficients
    A_use     = A_all(to_use,:);
    phase_use = input_phase(to_use);
    amp_use   = input_amplitude(to_use);
    
    %fit the good data
    model = NonLinearModel.fit(in, model_data, @SSD_fit_global, [input_period, phase_use, amp_use, A_use(:)'], 'Exclude', loc_same);%, 'Options', opts);
    
    %extract parameters from model
    Period          = model2array(model.Coefficients(1,1));
    PeriodError     = model2array(model.Coefficients(1,2));
    phases_good     =  model2array(model.Coefficients(2:number_used+1,1));
    phases_err_good = model2array(model.Coefficients(2:number_used+1,2));
    amps_good       = model2array(model.Coefficients(number_used+2:2*number_used+1,1));
    amps_err_good   = model2array(model.Coefficients(number_used+2:2*number_used+1,2));
    
    A_keep          = model2array(model.Coefficients(2*number_used+2:end,1));
    A_keep_err      = model2array(model.Coefficients(2*number_used+2:end,2));
    
    A_keep          = reshape(A_keep, size(A_use));
    A_keep_err      = reshape(A_keep_err, size(A_use));
    
    Resid_keep = model2array(model.Residuals(:,1));
    Fit_keep   = model2array(model.Fitted);
    
    for x = 1:numel(to_use)
        %points_used(x) = sum(whichbox(:)==to_use(x) & loc_same(:) ~= 1);
        %RChiSq(x) = ReducedChiSquared(model_data(whichbox==to_use(x)), Fit(whichbox==to_use(x)),loc_same(whichbox(:) == x));
        points_used(x) = sum(whichbox(:)==x & loc_same(:) ~= 1);
        RChiSq(x) = ReducedChiSquared(model_data(whichbox==x), Fit_keep(whichbox==x),loc_same(whichbox(:) == x));
    end
else
    Period = input_period;
    PeriodError     = NaN;
end


%generate output arrays of the correct size.
Phases          = NaN(1,number_boxes);
PhasesError     = NaN(1,number_boxes);
Amplitudes      = NaN(1,number_boxes);
AmplitudesError = NaN(1,number_boxes);
Ssd             = NaN(1,number_boxes);

% fit unused boxes and put data in correct places for good boxes.
for x = 1:number_boxes
    
    if ismember(x, to_use) %if boxes was used in period fit.
        %find location of the corresponding number
        loc = find(x == to_use);
        
        %copy numbers into the output array
        Phases(x)          = phases_good(loc);
        PhasesError(x)     = phases_err_good(loc);
        Amplitudes(x)      = amps_good(loc);
        AmplitudesError(x) = amps_err_good(loc);
        Ssd(x)             = nansum((Resid_keep(whichbox(:) == loc)).^2);
        
        A_all(x,:)         = A_keep(loc,:);
        A_allError(x,:)    = A_keep_err(loc,:);
        
        Residuals_temp(x,:)= Resid_keep(whichbox(:) == loc);
        fit_temp(x,:)      = Fit_keep(whichbox(:)==loc);
        
        n_points(x)        = points_used(loc);
        allChiSq(x)        = RChiSq(loc);
        
    else %if the box was not used in the period fit, fit the data.
        X = image_time;
        Y = pos_array(:,:,x);
        Z = SSD_array(:,:,x);
        
        X = X(:);
        Y = Y(:);
        Z = Z(:);
        A = A_all(x,:);
        
        in = [X Y];
        
        %find same image comparison data to exlude from the fitting.
        loc_same = (Z == 0);
        
        model = NonLinearModel.fit(in, Z, @SSD_fit_indv_phases, [input_phase(x), input_amplitude(x), A(:)'], 'Exclude', loc_same);%, 'Options', opts);
        
        Phases(x) =  model2array(model.Coefficients(1,1));
        PhasesError(x) = model2array(model.Coefficients(1,2));
        Amplitudes(x) = model2array(model.Coefficients(2,1));
        AmplitudesError(x) = model2array(model.Coefficients(2,2));
        Ssd(x) = model.SSE;
        
        A_all(x,:) = model2array(model.Coefficients(3:end,1));
        A_allError(x,:) = model2array(model.Coefficients(3:end,2));
        
        Fit = model2array(model.Fitted);
        
        Residuals_temp(x,:) = model2array(model.Residuals(:,1));
        fit_temp(x,:) = model2array(model.Fitted);
        
        n_points(x) = numel(Z) - sum(loc_same == 1);
        allChiSq(x) = ReducedChiSquared(Z, Fit, loc_same);
    end
end

%if required plot stuff.
if plots_on == 1
    
    surface_plot2(image_time, pos_array, SSD_array, fit_temp, Period)
    pause
%     for x = 1:number_boxes 
%         surface_plot(image_time, pos_array(:,:,x), SSD_array(:,:,x), fit_temp(x,:), Period)
%         if 1
%             pause
%         else
%             keyboard
%         end
%     end
end

%% vaidation of the fit (for data validation during debugging)
%Check that the residuals do not show large peak in fft.
if plots_on == 1%(0)
    disp('FFT of residuals for data validation.')  

    Z = reshape(Residuals_temp, size(SSD_array));
    Phases_fft( Z, image_time, 'poly', 'plotfigs');
    pause
end


%% combine data and prepare for output

%clean up phases by removing large values 
[Phases, Amplitudes] = Sinusoid_tidy(Phases, Amplitudes);

%Details for correcting reported lengths for the phase angle at the reference image. 
%determine the time of the reference image from analysis type.
if strcmpi(vars.analysis_type,'ref')
    ref_id = round(dat_size(2)/2);
elseif strcmpi(vars.analysis_type,'seq')
    ref_id = 1;
else
    error('The requested process type is not recognised')
end
RefDisp = Phases_SineFunc(image_time(1,ref_id), Period, Amplitudes, Phases);


%combine the boxes.
[PHASES, AMPLITUDES, reference_length, horizontal_position_pixels, PHASESerror, AMPLITUDESerror] = ...
    CombineBoxes(vars.length_type, 'phases', Phases, Amplitudes, data.boxX, data.boxY,...
                 'PhaseError', PhasesError, 'AmplitudeError', AmplitudesError, ...
                 'get_rid', vars.get_rid, 'symm_type', vars.symm_type,...
                 'RefOffsets', RefDisp);

%clean up phases by removing large values 
[PHASES, AMPLITUDES] = Sinusoid_tidy(PHASES, AMPLITUDES);

% get horizontal position of boxes (mid point of box)
horizontal_position_mm = horizontal_position_pixels * vars.scaling / 1000;

%removes variables from vars array because they are not used to
%cleans up the output files.
if strcmpi(vars.length_type, 'single') ~= 1 %if the boxes are pairs rename variable.
    number_pairs = length(PHASES);
    %     number_boxes = [];
end
if strcmpi(vars.length_type, 'pair') ~= 1 %if the boxes are not pairs remove the variables specified for pairs.
    vars = rmfield(vars, {'get_rid', 'symm_type'});
end


%% plot stuff
if 1
    figure(1)
    subplot(1,2,1)
    errorbar(PHASES, PHASESerror, 'x-')
    % errorbar(combined_phase, combined_phase_error)
    title('Phases')
    ylabel('Phase (radians)')
    xlabel('Radial position (boxes number)')
    
    subplot(1,2,2)
    errorbar(AMPLITUDES, AMPLITUDESerror, 'x-')
    % errorbar(combined_amplitude, combined_amplitude_error)
    ylabel('Amplitude (pixels)')
    xlabel('Radial position (boxes number)')
    title('Amplitude')
    
end

%% write output file

if output_save == 1
    
    disp('Writing the output files')
    
    % if the boxes are combined write two output files one of single valies
    % and other of combined values.
    
    %make structure of variables to write to output file.
    analysis_vars = vars;
    analysis_vars.polysurfcoeff = [m n];
    % output_vars.A = A_all;
    analysis_vars.fixed = fixed;
    
    %write single boxes to file if it has not already been done.
    clear output_values;
    single_fits.number_boxes        = {number_boxes};
    single_fits.horizontal_position = {mean(data.boxX, 2) * vars.scaling / 1000, 'mm'};
    single_fits.period              = {Period, 's'};
    if fixed == 0
        single_fits.period_err      = {PeriodError, 's'};
    end
    single_fits.phases              = {Phases, 'rad'};
    single_fits.phases_err          = {PhasesError, 'rad'};
    single_fits.amplitudes          = {Amplitudes, 'pixels'};
    single_fits.amplitudes_err      = {AmplitudesError, 'pixels'};
    single_fits.SSD                 = {Ssd, 'pixl^2'};
    single_fits.nos                 = {n_points, ''};
    single_fits.reference_position  = {mean(data.boxY, 2), 'pixels'};
    single_fits.useable             = {ismember(1:number_boxes, to_use), ''};
    single_fits.Chisq               = {allChiSq, ''};
    
    write_file_name = [out_file_name,'_single.txt'];
    
    analysis_vars.length_type = 'single';
    
    %write single foil fit file.
    WriteSineFit(write_file_name, single_fits,  analysis_vars, 'poly surf', A_all);
    
    %Write mat file with fit properties saved as matricies.
    save(Fits_file_name, 'A_all', 'A_allError', 'single_fits', 'analysis_vars', 'image_time', 'pos_array', 'SSD_array', 'image_time0');
    
    %     if x == 1
    if ~strcmpi(vars.length_type, 'single')
        
        %make structure of values to write to output file.
        combined_fits.number_pairs        = {number_pairs};
        
        combined_fits.horizontal_position = {horizontal_position_mm, 'mm'};
        combined_fits.period              = {Period, 's'};
        if fixed == 0
            combined_fits.period_err      = {PeriodError, 's'};
        end
        combined_fits.phases              = {PHASES, 'rad'};
        combined_fits.phases_err          = {PHASESerror, 'rad'};
        combined_fits.amplitudes          = {AMPLITUDES, 'pixels'};
        combined_fits.amplitudes_err      = {AMPLITUDESerror, 'pixels'};
        %         output_values.SSD                 = {Ssd, 'pixels.^2'};
        %         output_values.nos                 = {numel(SSD_array(:,:,1)), '#/foil'};
        if exist('reference_length', 'var')
            combined_fits.reference_length = {reference_length, 'pixels'};
        else
            combined_fits.reference_position = {reference_position, 'pixels'};
        end
        write_file_name = [out_file_name,'.txt'];
        
        WriteSineFit(write_file_name, combined_fits,  analysis_vars);
        
        save(Fits_file_name, 'combined_fits', '-append');
        %Combined_fits = output_values;
        
    end
end
% out = toc;
toc;

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


    function predict = SSD_fit_global(variables, predictors)
        
        num_boxes = predictors(1,4);
        M = predictors(2,4);
        N = predictors(3,4);
        t = predictors(:,1);
        displ = predictors(:,2);
        where = predictors(:,3);        
        
        per = variables(1);
        phas = variables(2:num_boxes+1);
        amp = variables(num_boxes+2:2*num_boxes+1);
        polycoef = variables(2*num_boxes+2:end);
        %polycoef = reshape(polycoef, polyvalnm_ncoef(M,N),num_boxes);
        polycoef = reshape(polycoef, num_boxes, polyvalnm_ncoef(N,M));
        
        for i = 1:num_boxes
%             predict(where==i,1) = SSD_model_primary(displ(where==i),t(where==i),per,phas(i),amp(i),polycoef(:,i)); 
            predict(where==i,1) = SSD_model_primary(displ(where==i),t(where==i),per,phas(i),amp(i),polycoef(i,:));        
        end
        
    end %SSD_fit_global


    function predict = SSD_fit_fixedT(variables, predictors)
        
        phas = variables(1);
        amp = variables(2);
        
        per = predictors(1,3);
        time = predictors(:,1);
        displ = predictors(:,2);
        nums_coef = predictors(2,3);
        polycoef = predictors(3:3+nums_coef-1,3);
        
        predict = SSD_model_primary(displ,time,per,phas,amp,polycoef);
        
    end %SSD_fit_indiv2

    function predict = SSD_fit_indv_phases(variables, predictors)
        per = Period;
        phas = variables(1);
        amp = variables(2);
        polycoef = variables(3:end);
%         polycoef = reshape(polycoef, polyvalnm_ncoef(n,m));
        
        t = predictors(:,1);
        displ = predictors(:,2);

        predict = SSD_model_primary(displ,t,per,phas,amp,polycoef);

    end %SSD_fit_indv_phases

    function [out] = SSD_model_primary(displacement, times, period, phase, amplitude, A)
        
        %displacement is the pixel displacement of the SSD value being calculated.
        % times - - time stamp of frame
        % period, phase and amplitude -- all of the driving wave
        % A -- polynomial surface coefficients. 
                
        %adjust x values based on the phase and ampltidue of the driving sine wave
        disp_min = Phases_SineFunc(times, period, amplitude, phase);
%         disp_min = amplitude .* cos(times./period*2*pi + phase);

        displacement = displacement + disp_min;
%        displacement = displacement - disp_min;
        %this function is a subtract not an addition because to collapse
        %the positions onto the background one has to remove the
        %contribution of the sine wave. The modified displacements are
        %those that would have been observed if there was no driving
        %wave present.

        %reshape polynomial coefficients
        coef = polyvalnm_coef2mat(A,n,m);
        
        %caclulate polynomial surface using adjusted x values.
        out = polyvalnm(coef,displacement,times);
        
    end %SSD_model_primary

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


%% VERSION info:
%  - v.2.5.1 - 9th April 2019
%   - Fixed Chisquared calculation and the distribution of calculated values which failed when not all boxes were used to fit period. 
%  - v.2.5 - 19th October 2018
%   - changed reference length so that now accounts for the displacment of the reference by the driving wave. 
%  - v.2.4 - 13 July 2018
%   - added 'no save' option
%  - v.2.3.1 - 5 July 2018
%   - fixed a bug in SSD_fit_global that was reordering the SSD surface
%   coefficients incorrectly. Now 
%  - v.2.3.1 - 5 June 2018
%   - changed components related to saving A_all which were passing the
%   wrong numbers around.
%  - v2.3 - 19 December 2017
%   - changed SSD_model_primary from 'displacement + min' to
%   'displacement-min' so that the fits correspond to those of Phases4_1.m
%  - v2.2 - 19 October 2017
%   - modified so that writes polynomial surface coefficients to *.mat file.
%  - v2.1 - 6 October 2017
%   - modified so that writes polynomial surface coefficients to the file
%  - v2.alpha - 3 July 2017
%   - Replaced the sine function with a call to Phases_SineFunc.m. this is
%   the same function that will be used by Phases.m and is now used because
%   the two functions had different formats in their sinusiod functions.
%   See Phases_SineFunc for details. alpha version is for testing.
%  - v1.8 - 26 June 2017
%   - Changed the input optiosns so that the period from the title can be
%   used expliciatlly -- this helps with batch processing.
%  - v 1.7.7 -- 23rd June 2017
%     - Reordered data preprocessing so that time stamp smoothing is after
%     cutting.
%  - v 1.7.6 -- 8th June 2016
%     - Added exclusions to the Chi Squared of the data that was exluded from the fitting.
%  - v 1.7.5 -- 6th June 2016
%     - Fixed the SSD and number of foils in the outputs so that the files can be read subequently.
%     - Added Chi Squared to outputs. 
%     - Minor changes to make sure the 'plot' option works without failing.
%  - v 1.7.4 -- 28th May
%     - Added 'out' option to file so that it will work with parfeval.
%     - Added SSD and number of foils to the outputs to be saved. 
%  - v 1.7.3 -- 20th May
%     - Fixed bug in image_time array setup. It now has to be correct
%     whichever way round the array is when it is read in. The script now
%     gives the same answers as v1.4.2 (which was a good version).
%  - v 1.7.2 -- 17th May
%     - Added options for setting polynomial order of surface to input
%     options.
%  - v 1.7.1 -- 17th May
%     - Changed to data read function to call ImageAnalysis_SSDio
%  - v 1.7 -- 7th May
%     - Edited the help text at the top of the file
%     - Added switch for plotting the figures if required.
%  - v 1.6 -- 9 March 2016
%     - Replace fourier transform section of code with a new stand alone
%     function called Phases_fft. This function will be compatible with the
%     fft function in Phases as well.
%     - Change output so that it always writes the single boxes data to a file.
%     - Removed data period validation to fft script. 
%     - reordered how the 'fixed' period is solved for.
%     - Moved data plotting into its own subfunction called surface_plot
%  - v 1.5 -- 26 Feb 2016
%     - Replace combining of boxes with a new stand alone function.
%  - v 1.4.2 -- 28 jan 2016
%     - bug fix for number of rows. Now finds new rows as long as there is
%     some overlap between the boxes.
%  - v1.4.1 -- 13 Jan 2016
%     - changed the say the outfile name is constructed to have a root and
%     a tail. This allows the Experiment_data_processing script to access
%     both parts of the name separately.
%  - v1.4 -- 12 Jan 2016
%     - Changed the input syntax to allow the period to be fixed to a value
%     and not solved for. 
%     - Made changes to allow a fixer period to be propagated through the
%     fitting routine.
%  - v1.3 -- 1 Jan '16
%     - Replaced difference and comdine sinusoid calls with new (properly tested) functions. 
%     - Replaced code tidying up of sine waves after fitting with a new
%     function called Sinsoid_tidy.
%  - v 1.2.3 -- 8th December 2015
%     - Bug fix. Changed order of lines associated with 'experiment_data_processing' so that empty data 
%     set does not crash the script.
%  - v1.2.2 -- 7 Dec 2015
%     - Changed how exclude SSD=0 data from fitting. Now use 'exclude' option in nonlinearmodel.
%  - v1.2.1 -- 7 Dec 2015
%     - Some tidying up of the script and removal of unused code.
%  - v1.2 -- 2 Dec 2015
%     - If the time stamps are granualar (ie. with 1 s precision) I now smooth them.
%     - added function which takes data processing unique to each
%     experiment out of this file and into a script called 'Experiment_data_processing'
%  - v1.1 -- ?? 
%     - Not sure of the change.
%     - possible changed the 'pair' case so that is not reports length
%     changes rather than whatever it was doing before that was wrong --
%     i.e. replaced CombineSineWaves with DifferenceWineWaves.
%  - v1.0 -- 18th October 2015
%     - Calculate the phases and amplitudes of the foil displacements from
%       the SSD used by ImageAnalysis/DisplacementAnalysis scripts.
%  - v1.1 -- 20th October
%     - Included the removal of no displacement comaprison (when SSD = 0)
%     - Changes to allow variable amount of displacements in analysis and
%       to cut the data size to match the search value.