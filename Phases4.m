% PHASES.
%   For calculating the phase of the standing wave for thermal diffusivity or 
%   anelasticity experiments. The script fits a sineusoid to the box
%   dispalcement data in '*_posiiotn_change.txt' seriies files returned by the 
%   ImageAnalysis.m script
%
%   syntax: Phases4(displaements_file_name.txt, options)
%     options: 'period', value  -- estimated value for period
%              'fixed'          -- force the period be that given by
%                        period. If fixed is set but no period set then fixed period
%                        is got from the file name.
%              'plot'           -- plot all the output and validation graphs.
%
%   The script is currently renamed because Phases v4+ uses Phases_SineFunc
%   as the sinusoid function which has different equation form to that used
%   in Phases v3.n and lower. See Phases_SineFunc for the deatils.
%
%   See also: ImageAnalysis, PhaseCheck, PhasesSSD, Phases_SineFunc

%   $ version: 4.aplha $ 4 July 2017 $
%   Simon Hunt 2008 - 2017
%       - for details of changes between versions see end of file.
%
%	- There was the begining of an option called 'combine' foils -- bit like pair but oposit. 
%  but this has been lost in all the rewriting.
%

function Phases4(file_name, varargin)
tic


%% settings and variables.

%check input file exists.
if ~exist(file_name, 'file')
   error('Phases:nofile', ['The input file, ', file_name, ', does not exist on the path. \n Check file name string and location.']) 
end

%set default options
fixed = 0;
input_period = NaN;
plots_on = 0;
leave = 0;

% prase varargin...
if length(varargin) == 1 && sum(isfinite(varargin{1})) == 1
    %this is a catch for the historical syntax. which expected just a
    %number for the period guess...
    warning('Phases:oldsyntax','Inputting only the a number for the period is old syntax')
    input_period = varargin{:};
else
    %parse new syntax.
    iarg = 1;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'period'
                input_period = varargin{iarg+1};
                iarg = iarg + 2;
            case 'fixed'
                fixed = 1;
                iarg = iarg + 1;
            case 'plot'
                plots_on = 1;
                iarg = iarg + 1;
                warning('the plotting bits have not been tested and are not complete')
            case 'exit'
                leave = 1;
                iarg = iarg + 1;
            otherwise
                error(['Unknown option: ' varargin{iarg}]);
        end
    end
end


% get experiment options and if not then make them.
expt = FileTitleInformation(file_name);
varfilename = char(strcat(expt(1), '.mat'));

%list of variables need to get/generate.
filevariables = {'length_type', 'get_rid', 'symm_type', 'bg_type', 'excessive', 'err_calc', 'scaling'};

%reads the experiment analysis options file. -- the file is called 'expt_name'.mat
if exist(varfilename, 'file')
    vars = load(varfilename, filevariables{:});
end

% if the variables are not read from 'expt_name'.mat run AnalysisOptions to generate them.
if sum(isfield(vars,filevariables)) ~= (length(filevariables))
    AnalysisOptions;
    vars = load(varfilename, filevariables{:});
end

%make filename for output file
[~, root, ~] = fileparts(file_name);
root_length = strfind(root,'_position_change');
root = root(1:root_length-1); %removes 'position_change' from the end of the input data file name.
tail = '_sine_fits_v4.txt';
out_file_name = [root,tail]; 

% minimisation options
minimise_options = optimset('MaxFunEvals',1E9, 'MaxIter',1E6 ,'TolX', 1e-3, 'Tolfun', 1e-3);

warning off all

%parse period given options
if strcmpi(input_period, 'title')==1
    input_period = expt{4};
    
elseif ~strcmpi(input_period, 'title') && ~isnumeric(input_period)
    error(['Error parsing inputs.\n The period is ''', input_period, ''' is not recognised.']);
    
elseif fixed == 1 && ~isfinite(input_period)
   error('Error parsing inputs.\n The period is fixed but there is no value to fix it to.') 
       
end

%% read data and initial processing

fprintf(1, '\n Phases of %s \n', file_name);

[foil_change_data, image_time, box_positions, ~, number_images] = ReadPositionChangeFile(file_name);

%moves the time stamps so that t=0 for the first image
image_time = image_time - image_time(1);

boxY = box_positions(1:2,:)';
boxX = box_positions(3:4,:)';

% Perform data pre-processing unique to experiment.
% the file 'Experiment_data_processing.m' needs to be in the current directory.
% It contains processing unique to the experiment e.g. for the Zn_? series of experiments the cutting of the begining and end of the data set is done here. 
if exist('Experiment_data_processing.m','file')
    preprocfunc = [mfilename,'_DataCut'];
    [foil_change_data, ~, image_time] = Experiment_data_processing(preprocfunc, foil_change_data, 0, image_time);

    %If too much has been cut exit otherwise recalcualte the dimension arrays, just incase they have changed.
    if numel(foil_change_data) == 0
        fprintf('No data left after data size reduction. Exiting \n\n')
        return
    else
        number_images = size(foil_change_data,1);
        %number_boxes = size(box_positions,1);
    end
else
    fprintf('Experiment_data_processing does not exist on path. No modifications to the data set.\n')
end


%% length change between pairs of boxes -- sorted by rows

%combine the boxes.
[length_change, image_time, reference_length, horizontal_position_pixels] = ...
    CombineBoxes(vars.length_type, 'displacements', foil_change_data, image_time, ...
                    boxX, boxY);

number_fits_todo = size(length_change,2);

if strcmpi(vars.length_type,'single')
    number_boxes = size(length_change,2);
    reference_position = reference_length;
    clear reference_length
else
    number_boxes = [];
    number_pairs = size(length_change,2);
end


% get horizontal position of boxes in mm 
horizontal_position_mm = horizontal_position_pixels * vars.scaling / 1000;
% starting_lengths_mm = reference_length * scaling / 1000;

% %correct for the 30degree offset in alignment of the foil
% %(NO LONGER NEEDED BUT LINE OF CODE REMAINS JUST INCASE).
% horizontal_position_mm = horizontal_position_mm / cos(30/180*pi);

%% initial period by fft

if isfinite(input_period)==0 && fixed == 0
    
    % Create initial guess for period by fft if one is not provided from varargin
    disp('FFT of data for rough period')
    
    %Call Phases_fft to calcualte the most likely period.
    arg_pass = {'old'};
    if plots_on == 1
        arg_pass = {arg_pass{:}, 'plotfigs'};
    end
    [ period bigger smaller input_phase ] = Phases_fft( length_change', image_time', arg_pass{:});
    
else
    %if not getting the range and phases from fft need to define them.    
    period = input_period;
    
    %define range of periods.
    bigger  = input_period*1.02;
    smaller = input_period*0.98;
    
    %define input phases.
    input_phase = ones(number_boxes,1);
end


%% refine period guess
% refine the period by minimising the range of the data as a function of phase.
% this is required because if the fft value is not close enough to the
% actual period the correct sine wave fitting does not work correctly.
%
% again if there is a period input to the calculation this step is missed out.

if exist('input_period','var') == 0   %FIX ME -- this section is now never called.

    spacing = (bigger-smaller)/100;
    periods = (smaller:spacing:bigger)';

    for x = 1 : length(periods)
        phase_array(:,x) = mod(image_time,periods(x));
    end

    bins = 7; %this number has to be prime -- to remove possibility of rational fractions of the number of frames per cycle (I think)

    for i = 1:number_fits_todo
        for x = 1 : length(periods)
            ranges = 0: periods(x)/bins :periods(x);
            %   plot(phase_array(:,x),detrended(:,pair),'r.')
            for y = 1: length(ranges)-1
                positions = find(phase_array(:,x) >= ranges(y) & phase_array(:,x) < ranges(y+1));
                if isempty(positions) == 1
                    spread(y,x) = NaN;
                else
                    values = detrended(positions,i);
                    spread(y,x) = max(values) - min(values);
                end
            end
            over(x) = sum(spread(:,x));
            over_all(i,x) = over(x);
        end
    end
    % for x = 1 : number_pairs
    %     min_range = find(over_all(x,:) == min(over_all(x,:)));
    %     minrange(x) = min_range(1);
    % end
    % consistent_min = mode(minrange);
    % period = periods(consistent_min);

    min_min_range = min(min(over_all));
    [x_min_pos,y_min_pos] = find(over_all == min_min_range);
    period = periods(y_min_pos(1));

%     % manually select the period from the sweep.      
%     % use incase period refining doesn't work
%     guess = 0;
%     guess2 = period;
%     while period ~= guess
%         close all;
%         plot(periods,over_all, '-', [guess2, guess2],[min(min(over_all)), max(max(over_all))], 'r');
%         pause;
%
%         guess = guess2;
%
%         prompt = {'Guess of the period'};
%         dlg_title = 'period guess';
%         num_lines = 1;
%         def = {num2str(guess2)};
%         answer = inputdlg(prompt,dlg_title,num_lines,def);
%
%         guess = str2num(answer{1,1});
%     %    period = guess;
%
%         if guess2 == guess
%             period = guess;
%             break
%         else
%             guess2 = guess;
%         end
%     end
%     period = guess
end


%% detrend the data
%removes background from data by deducting a moving average (of type specified) from the data.
switch vars.bg_type
    case 'moving average' %moving average
        %The length of the average is defined by the modal result from the fft and is approximately two periods long.
        average_over = 2*round(number_images/best_guess) + 1;
        if round(average_over/2) == average_over/2
            average_over = average_over + 1;
        end

        offset = average_over/2 + 0.5;
        image_time = image_time(average_over/2 + 0.5 : end - (average_over/2 + 0.5) + 1);

        for y = 1 : number_images - average_over + 1
            bg(y,:) = mean(length_change(y:y+average_over-1,:));
            detrended(y,:) = length_change(y + average_over/2-.5,:) -  bg(y,:);
        end
        
    case 'spline' %smothing spline
        %Smoothing takes place over a range defined by steps, which is the number of blocks the data
        %is divided into for the spline fitting.
        steps = round(number_images/100*3);
        for x = 1:number_fits_todo
            pp1 = splinefit(image_time,length_change(:,x),steps,3);
            spline_trend = ppval(pp1,image_time);
            detrended(:,x) = length_change(:,x) - spline_trend;

            if plots_on == 1
                bg(:,x) = spline_trend;
            end
        end

    case 'polynomial'
        for x = 1:number_fits_todo
            poly_average_coef = polyfit(image_time,length_change(:,x),2);
            detrended(:,x) = length_change(:,x) - polyval(poly_average_coef,image_time);
            
            if plots_on == 1
                bg(:,x) = polyval(poly_average_coef,image_time);
            end
        end

    otherwise
        error('Background removal method is not recognised, change something and try again')
end

%if requested plot the backgrounds.
if plots_on == 1
    
    if (1) %limit the number of plots that can be made.
        max_num_plot = 6;
        if number_fits_todo>max_num_plot
            disp('We are only plotting a subset of the data to save time')
            fit_plot = sort(unique(round(rand(1,max_num_plot)*number_fits_todo)));
        else
            fit_plot = 1:number_fits_todo;
        end
    else %plot all the fits! This could be a lot of plotting.
        fit_plot = 1:number_fits_todo;
    end
        
    %plot the selected data
    plot_fig = figure;
    for x = 1:length(fit_plot)
        
        plot_here = fit_plot(x);
        clf;
        
        subplot(2,1,1),
        plot(image_time,length_change(:,plot_here),'b',image_time, bg(:,plot_here), 'r'); 
        ylabel('displacement');
        xlabel('time');
        title('data and background');
        
        subplot(2,1,2)
        plot(image_time,detrended(:,plot_here),'b');
        ylabel('displacement');
        xlabel('time');
        title('detrended data');
        
        
        if strcmpi(vars.length_type, 'single') == 1
            leader = 'Box';
        else
             leader = 'Pair'; 
        end
        suptitle([leader,' ', num2str(plot_here)]);
        pause
    end
    
    if leave == 1
        return
    end
end


%% sine fit to the data.
 
% phase_time = mod(image_time,period);
% phase_rad = phase_time/period*2*pi;

%phase guess - by finding the maximum values in the detrended data.
for x = 1:number_fits_todo
    max_val_pos(x) = find(detrended(:,x) == max(detrended(:,x)),1,'first');
end
% phases = (phase_rad(max_val_pos)-pi/2)';
phases = (image_time(max_val_pos)/period*2*pi - pi/2)';

% finds RMS amplitudes of the data and defines an amplitude guess for use in the fitting.
amplitudes = [sqrt(mean(detrended.^2,1)) * sqrt(2); zeros(1,number_fits_todo)];

%fits phase and amplitude guesses to the data -- no period refinement
for x = 1 : number_fits_todo
    inputs = [phases(x) amplitudes(1,x) amplitudes(2,x)];
%     fit_all_data = fminsearch(@SineSQ_SSD,inputs,minimise_options,phase_rad,detrended(:,x),vars.excessive);
    fit_all_data = fminsearch(@SineSQ_SSD,inputs,minimise_options,image_time,detrended(:,x),period,vars.excessive);

    phases(x) = fit_all_data(1); %pastes returned parameters back into original array
    amplitudes(1,x) = fit_all_data(2);
    amplitudes(2,x) = fit_all_data(3);

%     SineSQ_SSD(fit_all_data,phase_rad,detrended(:,x)), pause %calculates SSD for fit
%     [SSD(x) ~] = SineSQ_SSD(fit_all_data,phase_rad,detrended(:,x),vars.excessive); %calculates SSD for fit
    [SSD(x) ~] = SineSQ_SSD(fit_all_data,image_time,detrended(:,x),period,vars.excessive); %calculates SSD for fit

%    initial_sin = amplitudes(1,x) * sin(phase_rad - phases(x)) + amplitudes(2,x)* cos(2*(phase_rad - phases(x)));
%    plot(phase_rad,detrended(:,x),'g.',phase_rad,initial_sin,'k.'); axis tight; xlabel('Phase (radians)'), ylabel('Offset (pixels)')
%    pause;
end

% make sure that all the phase values reported are less than +/- 2pi.
phases = mod(phases,2*pi);

%assumes that the data with the smallest SSD are the 'best'. Uses these to
%solve for the 'proper' period.
cutoff = 0.4;
ordered = sort(SSD);
if length(ordered) < 10
    max_allowed_SSD = max(SSD);
elseif length(ordered)*cutoff > 25
    max_allowed_SSD = ordered(25);
else
    max_allowed_SSD = ordered(round(length(ordered)*cutoff)-1);
end

y = 1;
for x = 1 : number_fits_todo
    if SSD(x) <= max_allowed_SSD
        useable(x) = 1;
        amp_not_bad(1,y) = amplitudes(1,x);
        amp_not_bad(2,y) = amplitudes(2,x);
        phase_not_bad(y) = phases(x);
        detrended_not_bad(:,y) = detrended(:,x);
        y = y + 1;
    else
        useable(x) = 0;
    end
end
number_pairs_not_bad = length(phase_not_bad);


%defines whether to use the cos^2 overtone part of the fitting or just to fit sine to
%the data.
ratio = amplitudes(2,:)./amplitudes(1,:);
ratio_bar = median(abs(ratio));
ratio_bar = 1; %makes sure that the sine_sq function part is not used.

if ratio_bar >= .1
    MODELb = @SineSSD;
%     MODELc = @sine;

    inputs = [period, phase_not_bad, amp_not_bad(1,:)];
    king_maker = fminsearch(@PartialDomination,inputs, minimise_options, image_time, detrended_not_bad, vars.excessive);

    period = king_maker(1);
    phase_not_bad = king_maker(2:number_pairs_not_bad+1);
    amp_not_bad(1,:) = king_maker(number_pairs_not_bad+2:2*number_pairs_not_bad+1);
    amp_not_bad(2,:) = zeros(size(amp_not_bad(2,:)));

    amplitudes(2,:) = zeros(size(amplitudes(2,:)));
else
    MODELb = @SineSQ_SSD;
%     MODELc = @sineSQ;

    inputs = [period, phase_not_bad, amp_not_bad(1,:) amp_not_bad(2,:)];
    king_maker = fminsearch(@GlobalDomination,inputs, minimise_options, image_time, detrended_not_bad, vars.excessive);

    period = king_maker(1);
    phase_not_bad = king_maker(2:number_pairs_not_bad+1);
    amp_not_bad(1,:) = king_maker(number_pairs_not_bad+2:2*number_pairs_not_bad+1);
    amp_not_bad(2,:) = king_maker(2*number_pairs_not_bad+2:end);
end

phase_time = mod(image_time,period)/period*2*pi; %phase time from 0 to 2pi
phase_for_fit = 0:(2*pi/100):2*pi;

y = 1;
for x = 1 : length(useable)
    if useable(x) == 0
        if ratio_bar >= .1
            inputs = [phases(x) amplitudes(1,x)];
        else
            inputs = [phases(x) amplitudes(1,x) amplitudes(2,x)];
        end
        results = fminsearch(MODELb,inputs,minimise_options,image_time,detrended(:,x),period,vars.excessive);
        phases(x) = results(1);
        amplitudes(1,x) = results(2);
        if ratio_bar < .1
            amplitudes(2,x) = results(3);
        end
    elseif useable(x) == 1
        phases(x) = phase_not_bad(y);
        amplitudes(1,x) = amp_not_bad(1,y);
        amplitudes(2,x) = amp_not_bad(2,y);
        y = y+1;
    else
        error('Something has gone horribly wrong')
    end
end

for x = 1:number_fits_todo
    if ratio_bar >= .1
        coefficients = [phases(x) amplitudes(1,x)];
%         allcalc = MODELc(phase_time, phases(x), amplitudes(1,x));
        allcalc = Phases_SineFunc(image_time, period, amplitudes(1,x), phases(x));
    else
        coefficients = [phases(x) amplitudes(1,x) amplitudes(2,x)];
%         allcalc = MODELc(phase_time, phases(x), [amplitudes(1,x) amplitudes(2,x)]);
        allcalc = Phases_SineFunc(image_time, period, amplitudes(1,x), phases(x), amplitudes(2,x));
    end
    [SSD(x) nos(x)] = MODELb(coefficients, image_time, detrended(:,x), period, vars.excessive);
    %FIX ME -- this line lloks wrong in that it is called hte wrong
    %function.

    miss_fit = detrended(:,x)-allcalc;
    in = abs(miss_fit) < std(miss_fit)*vars.excessive;
    Chisq(x) = ReducedChiSquared(detrended(in,x),allcalc(in));

    if strcmpi(vars.err_calc,'yes') == 1
        %error calculation
        % this is incomplete because it only works for the sine part of the
        % fiting function.
        out = errors(image_time,detrended(:,x), period, coefficients);
        phase_err(x) = out(1);
        amp_err(1,x) = out(2);
        if numel(out) >=3
            amp_err(2,x) = out(3);
        end
            
        %makes sure that the errors are not massive - which they can be if the data
        %is bad. If they are too large replaces huge number with 3 (i.e. still a comparatively large number).
        if phase_err(x) > 3
            phase_err(x) = 3;
        end
    else
        phase_err = NaN(size(SSD));
        amp_err = NaN(size(SSD));
    end


    %if required plot the final fit to each box/pair of boxes.
    if plots_on == 1
        
        %get the outsized residuals -- i.e. data not fitted to.
        %     outside = 4;
        %     mean_missfit = mean(miss_fit);
        %     std_missfit = std(miss_fit);
        %     outsized = find(abs(miss_fit) > std_missfit*outside);
        %     cut = std_missfit*outside;
        %     number_large = length(find(abs(miss_fit) > std_missfit*outside));
        
        if ismember(x,fit_plot) == 1
            gcf(plot_fig);
            clf;
            subplot(2,1,1); %cla
                plot(phase_time,detrended(:,x),'r.', phase_time,allcalc,'b.');
                xlabel('Phase (radians)');
                ylabel('delta length');
                title(strcat('pair',num2str(x), ' - global fit results'));
                xlim([0 2*pi]);
            subplot(2,1,2); %cla
                plot(phase_time, miss_fit,'g.',[0 2*pi],[0 0],'k-');
                xlabel('Phase (radians)');
                ylabel('delta length');
                title(strcat('pair ',num2str(x), ' - detrended data')); xlim([0 2*pi]);
                hold on;
                %             plot(phase_time(outsized), miss_fit(outsized),'k.');
                hold off;
            pause
        end
    end
end

%forces all the sinusiods to have positive amplitudes and 0<phases <2*pi
for x = 1 : number_fits_todo
    [phases(x), amplitudes(1,x)] = Sinusoid_tidy(phases(x), amplitudes(1,x));
end


%% write output file

%make structure of variables to write to output file.
output_vars = vars;
output_vars.fixed = fixed;

%make structure of values to write to output file.
if ~isempty(number_boxes)
    output_values.number_boxes = {number_boxes};
else
    output_values.number_pairs = {number_fits_todo};
end
output_values.horizontal_position = {horizontal_position_mm, 'mm'};
output_values.period              = {period, 's'};
if fixed == 0
%     output_values.period_err          = {PeriodError, 's'};
    warning('The period error does not exist as a value. Edit code to make it.')
end
output_values.phases              = {phases, 'rad'};
output_values.phases_err          = {phase_err, 'rad'};
output_values.amplitudes          = {amplitudes(1,:), 'pixels'};
output_values.amplitudes_err      = {amp_err(1,:), 'pixels'};
if size(amp_err,1) == 2
    output_values.amplitudes_sq       = {amplitudes(2,:), 'pixels'};
    output_values.amplitudes_sq_err   = {amp_err(2,:), 'pixels'};
    warning('WriteSineFit will not save the sine^2 values. It needs to be rewritten to incude them')
end
output_values.useable             = {useable, ''};
output_values.SSD                 = {SSD, 'pix^2'};
output_values.nos                 = {nos, ''};
output_values.Chisq               = {Chisq, ''};
if exist('reference_length', 'var')
    output_values.reference_length = {reference_length, 'pixels'};
else
    output_values.reference_position = {reference_position, 'pixels'};
end
write_file_name = out_file_name;

WriteSineFit(write_file_name, output_values,  output_vars);

warning on

toc
end %Phases

%% ###################### subfunctions #########################

% function Y = sine(X,phase,amplitude)
% if length(phase) == 1
%     Y = amplitude .* sin(X - phase);
% else
%     Y = amplitude(ones(1,length(X)),:) .* sin( X(:,ones(1,length(phase))) - phase(ones(1,length(X)),:)  );
% end
% 
% end %sine
% 
% function Y = sineSQ(X,phase,amplitude)
% Y = amplitude(1) .* sin(X - phase) + amplitude(2) .* cos(2*(X - phase));
% end %sineSQ

function [diff_sq nos] = SineSSD(coeff,X,Y, per, excess)
phase = coeff(1);
amplitude = coeff(2);

% Y_fun = sine(X,phase,amplitude);
Y_fun = Phases_SineFunc(X, per, amplitude(1), phase);

residuals = Y_fun - Y;

deviation = excess*std(residuals,0,1);
deviation = repmat(deviation,length(X),1);
discard = find(abs(residuals) > deviation);
residuals(discard) = 0;

nos = length(X) - length(discard);
diff_sq = sum(sum(residuals.^2));
end %SineSSD

function [diff_sq nos] = SineSQ_SSD(coeff,X,Y,per,excess)
phase = coeff(1);
amplitude = coeff(2:3);
% Y_fun = sineSQ(X,phase,amplitude);
Y_fun = Phases_SineFunc(X, per, amplitude(1), phase, amplitude(2));

residuals = Y_fun - Y;

deviation = excess*std(residuals,0,1);
deviation = repmat(deviation,length(X),1);
discard = find(abs(residuals) > deviation);
residuals(discard) = 0;
diff_sq = sum(sum(residuals.^2));

nos = length(X) - length(discard);
if amplitude(1)/amplitude(2) < 0
    diff_sq = diff_sq * 10e10;
end
end %sineSQ_SSD

function [diff_sq nos] = PartialDomination(coeff, X, Y, excess)
number_pairs = (length(coeff) - 1) / 2;
period = coeff(1);
phase = coeff(2:number_pairs+1);
amplitude(1,:) = coeff(number_pairs+2:2*number_pairs+1);
% 
% X_use = mod(X,period)/period*2*pi;
% Y_fun = sine(X_use, phase, amplitude);
Y_fun = Phases_SineFunc(X, period, amplitude, phase);

residuals = Y_fun - Y;

deviation = excess*std(residuals,0,1);
deviation = repmat(deviation,length(X),1);
discard = find(abs(residuals) > deviation);
residuals(discard) = [];

nos = length(X);
diff_sq = sum(sum(residuals.^2));
end %PartialDomination

function [diff_sq nos] = GlobalDomination(coeff, X ,Y, excess) %function name is David's suggestion
number_pairs = (length(coeff) - 1) / 3;
period = coeff(1);
phase = coeff(2:number_pairs+1);
amplitude = coeff(number_pairs+2:number_pairs*2+1);
amp_sq = coeff(number_pairs*2+2:length(coeff));
number_times = length(X);

X_use = mod(X,period);
X_use = X_use/period*2*pi;


Y_fun = Phases_SineFunc(X, period, amplitude, phase, amp_sq);
% Y_fun = amplitude(ones(1,number_times),:) .* sin(X_use(:,ones(number_pairs,1)) - phase(ones(1,number_times),:))...
%      + amp_sq(ones(1,number_times),:) .* cos(2*(X_use(:,ones(number_pairs,1)) - phase(ones(1,number_times),:)));


% for x = 1 : number_pairs
%     Y_fun(:,x) = amplitude(x) * sin(X_use - phase(x)) + amp_sq(x) * cos(2*(X_use - phase(x)));
%   %  plot(X_use,Y_fun(:,x),'k.',X_use,Y(:,x),'g.')
%   %  pause
% end


residuals = Y_fun - Y;
% mean(residuals,1);
deviation = excess*std(residuals,0,1);
deviation = repmat(deviation,length(X_use),1);
discard = (abs(residuals) > deviation);
residuals(discard) = [];

nos = length(X);
diff_sq = sum(sum(residuals.^2));

if min(amplitude./amp_sq) < 0
    diff_sq = diff_sq * 10e10;
end
end %GlobalDomination

function se = errors(X_phase,y,per, coefficients)
% taken from err.m... need to annotate this part properly!

dof = length(X_phase) - length(coefficients);

% standard deviation of the residuals
sdr = sqrt(SineSSD(coefficients,X_phase,y,per,inf)/dof);

% jacobian matrix
% model_to_pass = @(coefficients) sine(X_phase, coefficients(1), coefficients(2));
model_to_pass = @(coefficients) Phases_SineFunc(X_phase, per, coefficients(1), coefficients(2));
J = jacobianest(model_to_pass,coefficients);

% I'll be lazy here, and use inv. Please, no flames,
% if you want a better approach, look in my tips and
% tricks doc.
Sigma = sdr^2*inv(J'*J);

% Parameter standard errors
se = sqrt(diag(Sigma))';
end %errors


function RedChi2 = ReducedChiSquared(obs,model)
%copied from chi2test on mathworks website.

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


%% Known issues
% -- The detrending by spline and presumably moving average does not work
% with data with only a few (about 10) cycles in it.

%% Versioning
% v3.2 - 26 June 2017
%   - Changed the input optiosns so that the period from the title can be
%   used expliciatlly -- this helps with batch processing.
% v3.1 - 23 June 2017
%   - Added data cutting from 'Experiment_data_processing' (if present) and
%   reordered the data reading the processing to match that of PhasesSSD.m
% v3.0 - 7th May 2016
%   - Improved help at top of file
%   - Changed the way the varialbe inputs are processed and added options. So can now fix the period to a value which is input.
%   - Moved fft for period to a separate script called phases_fft.
%   - Added a switch to plot the fits if requested.
%   - replaced the attempt to clean up the outputs with calls to Sinusois_Tidy.m
% v2.3 - 12 August 2014
%   - Changed the reading of the input files so it now uses script ReadPositionChangeFile script
%   - Changed the writing of the output files so WriteSineFit is now used. 
%   - Both these new scripts contain the functionality previously contained within this script.
% v2.2 - 26 April 2013
%   - Added functionality so can measure phases in data with single boxes in each row (i.e. mechanical deformation experiments).
% v2.1.2 - 25 Nov 2011
%   - Fixed error in fourier transform rourine which only seemed to be important when more boxes than images
% v2.1.1 - nov 2011
%   - changed number of rows to skip for image_times - following from fix in ImageAnalysis.m
% v2.1, December 2010
%   - Moved the analysis options to the file 'experiment_name'.mat and if they don't exist generates
%   them in AnalysisOptions.m
%   - Now reads m file verions number from header rather than separate variable, moved sine_fits
%   header generation to be with the rest of the out file generation.
% November 2010: (called Phases v2.0 -- version 1 was called PhaseLag)
%   - Added version numbers
%   - Added header to output file which lists the options used in the analysis - useful when reviewing the processed data.
%   - RENAMED FILE from PhaseLag (which it doesn't find) to PHASES (which it does)
%   - Now two sorts of background removal spline and moving average
%   - Removed the outliers (>n sigma) from displacement data when calculating the phase and amplitudes.
%   - Options to discard different foils in case of odd number thereof.
%   - Can now use individual foils as 'length changes'
% 21st June:
%   - Make sure that in the outfile amplitudes>0 and range(phases) < pi.
%     This aids automatic data processing later.
%   - Limits the size of the errors to 3 radians (a large, effectivly random,
%     number).
% 11th June 2010:
%	- incorporated various updates to the code (many of which I don't remeber) but it does
%	  include:
%		>	tidying up the solving routine and rationalising the number of subfunctions for sine solving
%		>	incorporating an option to find the phases between asymmettic boxes (used for CaIrO3 expt)
%		>	Removed the need for their to be an even number or foils in the original data it now assumes that
%				the largest gap between foils brackets the middle of the sample.
%		>	An error calculating routine has been added (at some stage) but it is very slow.
%		>	the rest of the code has been tidied up and cleaned of a number of spare and repeated variables.
%	*** The phase error calculation only works for the sine, not sine squared function
%           as well as being faily slow
% 13th Oct 2009:
%   - changed input so script can take input phase from elsewhere (usually this will be PhaseCheck.m).
% 3rd Oct 2009:
%   - added routine to make sure that there are an even number of rows fed to the fitting routine.