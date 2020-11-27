function [ detrended, BackGround, varargout ] = Phases_Detrend(det_type, length_change,  image_time, varargin)
%Phases_detrended. Detrend displaments for Phases Cancluaiton 
%
%   Syntax:
%     [detrended, background] = Phases_detrend(bg_type, length_changes, TimeStamps, [additional args])
%         bg_type:       background type to remove
%         length_changes:array of length changes
%         TimeStamps:    array of time stamps
%
%         recognised bacgfound types:
%         'moving average' - requires period as [additional arg]
%         'polynomial'     - optional argument for order of polynomial as [additional arg]
%         'spline'         - optional argument for number of segments as [additional arg]
% 
%   See also: ImageAnalysis, PhaseCheck, Phases
   
% Simon Hunt. June 2018.
%  Version 1




%% parse input
num_inputs = nargin;
number_images = length(image_time);
number_fits_todo = size(length_change,2);

plots_on = 0;
if ismember(varargin, 'plot')
    plots_on = 1;
    varargin(cellfun(@(x) strcmpi(x, 'plot'), varargin)) = [];
    num_inputs = num_inputs-1;
end

switch det_type
    case 'moving average'
        best_guess = varargin{1};
        
    case 'polynomial'
        if num_inputs > 3
            poly_order = varargin{1};
        else
            poly_order = 5;
        end
    case 'spline'
        if num_inputs > 3
            steps = varargin{1};
        else
            steps = round(number_images/100*3);
        end
    otherwise
        error('Background removal method is not recognised, change something and try again')
        
end

%% detrend the data
%removes background from data by deducting a moving average (of type specified) from the data.
switch det_type
    case 'moving average' %moving average
        %The length of the average is defined by the modal result from the fft and is approximately two periods long.
        average_over = 2*round(number_images/best_guess) + 1;
        if round(average_over/2) == average_over/2
            average_over = average_over + 1;
        end

        offset = average_over/2 + 0.5;
        image_time = image_time(average_over/2 + 0.5 : end - (average_over/2 + 0.5) + 1);

        for y = 1 : number_images - average_over + 1
            BackGround(y,:) = mean(length_change(y:y+average_over-1,:));
            detrended(y,:) = length_change(y + average_over/2-.5,:) -  BackGround(y,:);
        end
        
    case 'spline' %smothing spline
        %Smoothing takes place over a range defined by steps, which is the number of blocks the data
        %is divided into for the spline fitting.
        for x = 1:number_fits_todo
            pp1 = splinefit(image_time,length_change(:,x),steps,3);
            spline_trend = ppval(pp1,image_time);
            detrended(:,x) = length_change(:,x) - spline_trend;

            BackGround(:,x) = spline_trend;
        end

    case 'polynomial'
        for x = 1:number_fits_todo
            poly_average_coef = polyfit(image_time,length_change(:,x),poly_order);
            detrended(:,x) = length_change(:,x) - polyval(poly_average_coef,image_time);
            
            BackGround(:,x) = polyval(poly_average_coef,image_time);
        end
end

%if requested plot the backgrounds.
if plots_on == 1
    
    if (1) %limit the number of plots that can be made.
        max_num_plot = 6;
        if number_fits_todo>max_num_plot
            % radonmly get a subset of the boxes to plot
            disp('We are only plotting a subset of the data to save time')
            fit_plot = sort(unique(round(rand(1,max_num_plot-1)*number_fits_todo))+1);
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
        plot(image_time,length_change(:,plot_here),'b',image_time, BackGround(:,plot_here), 'r'); 
        ylabel('displacement');
        xlabel('time');
        title('data and background');
        
        subplot(2,1,2)
        plot(image_time,detrended(:,plot_here),'b', [min(image_time) max(image_time)], [0 0], 'r');
        ylabel('displacement');
        xlabel('time');
        title('detrended data');
        
        %length_type = evalin('caller', 'vars.length_type');
        %if strcmpi(length_type, 'single') == 1
            leader = 'Box';
        %else
             leader = 'Pair'; 
        %end
        suptitle([leader,' ', num2str(plot_here)]);
        pause
    end
    
end
