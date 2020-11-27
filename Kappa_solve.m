% Solves for thermal diffusivity.
% Simon Hunt 2009 - 2014
% Takes the output from Phase_lag.m and solves for the thermal
% diffusivity after the method of Dobson et al. 2008 (Min mag)
%
%   $ version: 1.8 $ 13th August 2014 $
%     - for details of changes between versions see end of file.

function Kappa_solve(varargin)

number_files = nargin;
file_names = char(varargin);

%% read and organise input data. 

% get temp, load and experiment name.
for x = 1 : number_files
    file = file_names(x,:);
    expt_info = FileTitleInformation(file);
    expt_name = expt_info{1};
    force(x) = expt_info{2};
    temp(x) = expt_info{3};
    period(x) = expt_info{4};
    run_num(x) = expt_info{5};
end

% sorts the files by period - in descending order (so the data is in the
% right order for later in the analysis
period_list(:,1) = period;
period_list(:,2) = 1:number_files;
period_sort = sortrows(period_list);
period_sort = flipud(period_sort);
order = period_sort(:,2);
file_names(1:number_files,:) = file_names(order,:);
filenames = varargin(order);

%read the data from the files.
qc = ReadSineFit('open',filenames{:});
if qc == 1
    %reads list of headers
        titles = ReadSineFit('headers')
    %reads radii from all_numbers
        all_radii_mm = ReadSineFit('position');%(:,:) = all_numbers(:,2,:);
    %reads periods from all_numbers    
        all_periods = ReadSineFit('Period');%(:,:) = all_numbers(:,3,:);
    %reads phases from all_numbers        
        all_phases = ReadSineFit('Phases');%(:,:) = all_numbers(:,4,:);
        all_phases_err = ones(size(all_phases))*1e-5;%ReadSineFit('Phases error');%(:,:) = all_numbers(:,4,:);
        if isnan(all_phases_err) == 1
            all_phases_err = ones(size(all_phases)) * 0.05;
            warndlg('There are no errors for the calculated phases. The error on each phase has been assumed to be 1');
        end
    %reads amplidues from all_numbers    
        all_amplitudes = ReadSineFit('amplitudes');%(:,:) = all_numbers(:,5,:);
    %reads lengths from all_numbers    
        all_lengths = ReadSineFit('lengths');%(:,:) = all_numbers(:,9,:);
    %reads SSD from all_numbers    
%         all_SSD = ReadSineFit('SSD');%(:,:) = all_numbers(:,8,:);    

    % makes weights from phase errors
        all_weights = 1./(all_phases_err.^2);
    
end
qc = ReadSineFit('close');

%flips the phase data. This has to be done to get the phase lag in the
%right orientation.
all_phases = -all_phases;

%reads the experiment analysis options file. -- the file is called 'expt_name'.mat
if exist(strcat(expt_name,'.mat'),'file') == 2
    load(strcat(expt_name,'.mat'), 'R_solve', 'R__nought','template');
end
if exist('R_solve','var') == 0
    R_solve = 1; %set to 1 to solve for r0, otherwise set to 0 to keep fixed value
end

%checks the data to assertain if there is more than one pair of foils (i.e.
%more than one set of phase lags) in the data read in.
diff = all_radii_mm(1,2:end) - all_radii_mm(1,1:end-1);
num_rows = length(find(diff < 0)) + 1;

ends = find(diff(1,:) < 0);
ends = [0 ends length(all_radii_mm)];

%checks to see if there are any negative amplitues - if there are it corrects them.
neg_amps = find(all_amplitudes < 0);
if isempty(neg_amps) == 0
    all_amplitudes(neg_amps) = -all_amplitudes(neg_amps);
    all_phases(neg_amps) = all_phases(neg_amps) + pi;
    %need to correct phase squared as well if it ever makes an apperance.
end
    
%sets minimum radius to zero
min_rad = min(all_radii_mm,[],2);
all_radii_mm = all_radii_mm - repmat(min_rad,1,size(all_radii_mm,2));

% other variables that need to be set.
spacing = 0.1;
colours = get(0,'defaultAxesColorOrder'); %this gets the default colour scheme for figures. It does the same as colormap(lines) but without creating the figure

%% Build Window
if isempty(findobj('Type','figure','Tag','KappaWindow'))
    %define sizes needed for display panels

 
    %define column sizes for panels.
    gap1 = 30;
    gap2 = 10;
    column = 60;
    height = 15;
    border = 3;
    
    column1 = border;
    column2 = column1 + column + gap2;
    column3 = column2 + column + gap2;
    pan_width = 2*column + gap2 + 4*border;
    
    %extra size for the display panel
    extra_width = pan_width; %pixels
    top_margin = 0; %pixels    
        
    %figure creation
    fig = figure('Name','Kappa Solve','NumberTitle','off', 'Tag', 'KappaWindow','Visible','off');
    dims = get(fig,'Position');  
    set(fig, 'Position',[dims(1) dims(2) dims(3)+extra_width dims(4)+top_margin]); %set figure with extra margins
    movegui(fig,'center'); % Move the figure to the center of the screen.    
    
    % Axes creation
    axes_def = get(0,'DefaultAxesPosition');
    axes_pos = [axes_def(1)*dims(3),axes_def(2)*dims(4),axes_def(3)*dims(3),axes_def(4)*dims(4)];
    kap_plot = axes('Units','Pixels','Position',axes_pos,'Tag','kappa_rad_plot');

    % axes labels
    set(get(kap_plot,'XLabel'),'String','Offset (mm)');
    set(get(kap_plot,'YLabel'),'String','Phase lag (rad)');
    box on
    
    %define last corner of graph.
    end_fig = [axes_pos(1)+axes_pos(3) axes_pos(2)+axes_pos(4)];
    
    %Experimental Information Panel.
    pan_info_height = 4*height + gap2 + 2*border;
    row = [height*3 height*2 height 0] + border;
    info_bl = [end_fig(1)+gap1 end_fig(2)-pan_info_height]; %position of the bottom-left corner of the experimental information panel
        %make panel
    panel_info = uipanel('Title', 'Experiment details', 'Units','Pixels', 'Position', [info_bl(1) info_bl(2) pan_width pan_info_height], 'Parent', fig);
        %panel components
    title{1} = uicontrol('Style','text','String','Name',       'Position',[column1,row(1),column,15],'FontWeight','bold','Parent',panel_info);
    title{2} = uicontrol('Style','text','String','Load (ton)', 'Position',[column1,row(2),column,15],'FontWeight','bold','Parent',panel_info);
    title{3} = uicontrol('Style','text','String','Temp (°C)',  'Position',[column1,row(3),column,15],'FontWeight','bold','Parent',panel_info);
    title{4} = uicontrol('Style','text','String','Foil pair #','Position',[column1,row(4),column,15],'FontWeight','bold','Parent',panel_info);
    info{1} = uicontrol('Style','text','Position',[column2,row(1),column,15],'Parent',panel_info,'Tag','info-1');
    info{2} = uicontrol('Style','text','Position',[column2,row(2),column,15],'Parent',panel_info,'Tag','info-2');
    info{3} = uicontrol('Style','text','Position',[column2,row(3),column,15],'Parent',panel_info,'Tag','info-3');
    info{4} = uicontrol('Style','text','Position',[column2,row(4),column,15],'Parent',panel_info,'Tag','info-4');

    % Kappa solution panel setup
    pan_height = (number_files+4+1)* height + border*2;
    row = (((number_files+4)-1:-1:0) * height) + border;
    panel_bl = [end_fig(1)+gap1 end_fig(2)-pan_height-gap2-pan_info_height]; %position of the bottom-left corner of the panel
        %make panel
    panel1 = uipanel('Title', 'Fitted values', 'Units','Pixels', 'Position', [panel_bl(1) panel_bl(2) pan_width pan_height], 'Parent', fig);
    title{5} = uicontrol('Style','text','String','r0 (mm):',      'Position',[column1,row(1),column,15],'FontWeight','bold','Parent',panel1);
    title{6} = uicontrol('Style','text','String','period (s)',    'Position',[column1,row(3),column,15],'FontWeight','bold','Parent',panel1);
    title{7} = uicontrol('Style','text','String','kappa (mm²s-1)','Position',[column2,row(3),column,15],'FontWeight','bold','Parent',panel1);
%     title{8} = uicontrol('Style','text','String','\Chi^2','Position',[column3,row(3),column,15],'FontWeight','bold','Parent',panel1);
    radius_0_disp = uicontrol('Style','text','Position',[column2,row(1),column,15],'Parent',panel1,'Tag','rad_disp');
    period_display{1} = uicontrol('Style','text','String','all','Position',[column1,row(4),column,15],'Parent',panel1);
    kappa_display{1} = uicontrol('Style','text','Position',[column2,row(4),column,15],'Parent',panel1,'Tag',['kap_disp-',num2str(1)]);
    for x = 1 : number_files
        row_num = 4+number_files-x+1;
        period_display{x+1} = uicontrol('Style','text','Position',[column1,row(row_num),column,15],'ForegroundColor',colours(x,:),'Parent',panel1,'tag',['per_disp-',num2str(x+1)]);
        kappa_display{x+1} =  uicontrol('Style','text','Position',[column2,row(row_num),column,15],'ForegroundColor',colours(x,:),'Parent',panel1,'tag',['kap_disp-',num2str(x+1)]);%,'HorizontalAlignment','right');
    end
    
    % % position and contents of instructions text
    % instruct = 'Left mouse button to remove data; right button to add; ''s'' to change spacing; ''z'' to change zoom; ''p'': add 2pi; ''m'': remove 2pi';
    % instruct_text = uicontrol('Style','text','String',instruct,'Position',[10,400,600,15],'BackgroundColor',[.8 .8 .8]);

    % Change units to normalized so components resize automatically.
    % set([kap_plot,panel1,panel_info,instruct_text],'Units','normalized')
    set([kap_plot,panel1,panel_info],'Units','normalized')
%     set(title,'Units','normalized')

    set(fig,'Visible','on'); %make figure visible

else
    %if window exists then use it
    % the following lines define the labels for the changing parts
    fig = findobj('Type','figure','Tag','KappaWindow');
    kap_plot = findobj('Type', 'axes', 'Tag', 'kappa_rad_plot');
    axes_pos = get(kap_plot, 'Position');
    
    info{1} = findobj('Tag', 'info-1');
    info{2} = findobj('Tag', 'info-2');
    info{3} = findobj('Tag', 'info-3');
    info{4} = findobj('Tag', 'info-4');
    radius_0_disp = findobj('Tag','rad_disp');
    kappa_display{1} = findobj('Tag',['kap_disp-',num2str(1)]);
    for x = 1:number_files
        period_display{x+1} = findobj('Tag',['per_disp-',num2str(x+1)]);
        kappa_display{x+1} =  findobj('tag',['kap_disp-',num2str(x+1)]);
    end
end


%% solve for kappa
for count = 1 : num_rows
    row_to_use = num_rows - count + 1; %which set of phases to use in the data
    which_set = num2str(count); %so can display in data bar on right of figure
   
    %selects data to use - depends in row to use. 
    radii_mm = all_radii_mm(:,ends(row_to_use)+1:ends(row_to_use+1));
    periods = all_periods(:,ends(row_to_use)+1:ends(row_to_use+1));
    phases = all_phases(:,ends(row_to_use)+1:ends(row_to_use+1));
    phase_err = all_phases_err(:,ends(row_to_use)+1:ends(row_to_use+1));
    lengths = all_lengths(:,ends(row_to_use)+1:ends(row_to_use+1));
%     SSD = all_SSD(:,ends(row_to_use)+1:ends(row_to_use+1));
    weights = all_weights(:,ends(row_to_use)+1:ends(row_to_use+1));
    
    %adjust the phase data... to remove any phase wrapping problems.     
    phases = -phases; %+ pi;
    phases = rem(phases,2*pi);
    range_all_phases = max(phases,[],2)-min(phases,[],2);

    %define a few other parameters
    number_pairs = length(radii_mm);
    omegas = (2*pi)./periods;
    phi0 = min(phases,[],2);
    filter = ones(size(phases));

    % checks for tempate and if not discards underconstrained values
    if exist('template','var') == 0
        %discards data with large very large SSD values
%         acceptable = 2;
%         mean_SSD = mean(SSD,2);
%         mean_SSD = repmat(mean_SSD,1,number_pairs);
%         normalised_SSD = SSD ./ mean_SSD;
%         large = find(normalised_SSD > acceptable);
%         filter(large) = 0;
        
%         %discards the end 'n' radii.
%         discard = 1;
%         filter(:,1:discard) = 0;
%         filter(:,size(filter,2)-discard+1:size(filter,2)) = 0;

        %discards data with large standard errors
        too_big = 0.2;
        large = find(phase_err > too_big);
        filter(large) = 0;
    else
        filter = template;   
    end

    % looks for initial r_0 guess if no value is set
    if R_solve == 1 || exist('R__nought','var') == 0
        for i = 1 : number_files
            try
                radii_for_guess = radii_mm(i,:);
                phases_for_guess = phases(i,:)- phi0(i);
                filter_for_guess = filter(i,:);
                
                guess = polyfit(radii_for_guess(filter_for_guess==1),phases_for_guess(filter_for_guess==1),2);
                estimate(i,:) = polyval(guess,radii_for_guess);
                
                min_est = polyder(guess);
                r0_guess(i) = -min_est(2)/min_est(1);
            catch
                r0_guess(i) = mean(radii_for_guess);
            end
            
        end
    else
        r0_guess = R__nought;
    end
    r0_mm = mean(r0_guess);

    %other initial values needed.
    kappa = 1;
    kappa_guesses = repmat(kappa,1,number_files);

    howmany = 1:number_files;
    period = howmany(ones(1,number_pairs),:)';
    phi_nought = phi0';

    limits_r = [min(min(radii_mm))-.05 max(max(radii_mm))+.05];
    limits_k = [0 100];
   
    if exist('y_scaling','var') == 1
        clear y_scaling
    end

    done = 0;
    while done == 0

        if kappa > .95 * limits_k(2);
            kappa = 1;
            phi_nought = phi_nought + .1;
        end
        
        %fits the all the data to one kappa value
        which_period = period(filter==1);
        good_omega = omegas(filter==1);
        phases_use = phases(filter==1);
        radii_use = radii_mm(filter==1);
        weights_use = weights(filter==1);
        
        if R_solve == 1 %if the r0 is to be solved for as well
            coeffs = [kappa r0_mm phi_nought];
            
%             bestcoeffs4kappa=fminsearch(@kappasearch,coeffs,optimset('FunValCheck','on','MaxFunEvals',1E9,'MaxIter',1E6),...
%                 radii_use,phases_use,good_omega,which_period,limits_r,limits_k);
            bestcoeffs4kappa=fminsearch(@kappasearchWEIGHTED,coeffs,optimset('FunValCheck','on','MaxFunEvals',1E9,'MaxIter',1E6),...
                radii_use,phases_use,good_omega,which_period,limits_r,limits_k,weights_use);
            
            kappa = abs(bestcoeffs4kappa(1));
            r0_mm = bestcoeffs4kappa(2);
            phi_nought = bestcoeffs4kappa(3:end);
            
            radii_calc = abs((radii_mm - r0_mm));
            phi_all = CalcPhase(radii_calc,omegas,kappa);
        else % if r0 is fixed.
            radii_use = abs(radii_use - r0_mm);
            inputs = [R_solve kappa phi_nought];
            estimates  = fminsearch(@KappaFind,inputs,optimset('FunValCheck','on','MaxFunEvals',1E9,'MaxIter',1E6),...
                radii_use,phases_use,good_omega,filter,limits_k,weights_use);
            
            kappa = abs(estimates(2));
            phi_nought = estimates(3:end);
            
            radii_calc = abs((radii_mm - r0_mm));
            phi_all = CalcPhase(radii_calc,omegas,kappa);
        end
        
        calc_lag_step = (radii_mm(1,2) - radii_mm(1,1))/50;
        calc_lag_r = radii_mm(1):calc_lag_step:radii_mm(end);
        calc_lag_ra = abs((calc_lag_r - r0_mm));
        calc_lag_ra = calc_lag_ra(ones(number_files,1),:);
        omegas_calc = omegas(:,ones(1, length(calc_lag_ra)));
        calc_lag_phi = CalcPhase(calc_lag_ra,omegas_calc,kappa);
                
        adjustment = phi_nought';
        adjustment = adjustment(:,ones(1,number_pairs));
        
        %finds the most significant data point in the combined sets.
        significance = abs(phases - phi_all)./sqrt(weights);
        significance = significance .* filter; % remove preremoved data from find
        [C,I] = max(significance);
        [C2] = max(C);
        most_significant = [find(C2==C), I(C2==C)];
        
        %calc Reduced Chi^2
        corr = phases(filter==1)-adjustment(filter==1);
        ChiSq_all = ReducedChiSquared(corr,phi_all(filter == 1),weights(filter == 1));
        fprintf(1, 'reduced Chi squared = %6.2f \n', ChiSq_all)
        
        %fits the data to the bessel function for each data set individually
        for j = 1 : number_files

            filter_for_guess = filter(j,:);

            if sum(filter_for_guess) ~= 0
                radii_for_guess = radii_mm(j,:);
                radii_for_guess = abs(radii_for_guess(filter_for_guess==1) - r0_mm);

                phases_for_guess = phases(j,:);
                phases_for_guess = phases_for_guess(filter_for_guess==1);

                omega = omegas(j,:);
                omega = omega(filter_for_guess==1);
                
                weights_guess = weights(j,:);
                weights_guess = weights_guess(filter_for_guess==1);
                

                inputs = [kappa_guesses(j) phi0(j)];
                try
                    estimates = fminsearch(@kappaonly,inputs,optimset('FunValCheck','on','MaxFunEvals',1E9,'MaxIter',1E6),...
                        radii_for_guess,phases_for_guess,omega,limits_k,weights_guess);

                
%                   CalcPhase(radii_for_guess, omega, kappa_guesses(j));%, pause
                
                    kappa_guesses(j) = abs(estimates(1));
                    phi0(j) = estimates(2);
                catch
                    warning(['''Kappaonly'' returned NaN: The offening period is ', num2str(2*pi/omega(1))])
                end
            end
        end
        phi = CalcPhase(radii_calc,omegas,(kappa_guesses(ones(1,number_pairs),:))');

        %plot the data and the fits.
        adjustment = phi_nought';
        adjustment = adjustment(:,ones(1,number_pairs));

        offsets_set = (0:spacing:(number_files-1)*spacing)';
        offsets = offsets_set(:,ones(1,number_pairs));

        Phases_to_plot = phases - adjustment + offsets;

        phi = phi + offsets + phi0(:,ones(1,number_pairs)) - adjustment;
        phi_all = phi_all + offsets;

%         offsets = offsets_set(:,ones 1);
        calc_lag_phi = calc_lag_phi + offsets_set(:,ones(1,length(calc_lag_r)));

        
        figure(fig);  %makes main figure active figure.
        set(kap_plot, 'NextPlot', 'add'); %-- hold figure without adding text to command window  
        cla(kap_plot);  %clear axes.
        
        if exist('y_scaling','var') == 1
            set(kap_plot,'Ylim',y_scaling);
        else
            ylim(kap_plot,'auto');
        end
        set(kap_plot,'Xlim',limits_r);     
            
        for k = 1 : number_files
            filter_plot = filter(k,:);
            phases_plot = Phases_to_plot(k,:);
            radii_plot = radii_mm(k,:);
            err_plot = phase_err(k,:);

            %these plots should probably use scatter rather than plot to do the plotting.
            plot(kap_plot,radii_plot(filter_plot==1)', phases_plot(filter_plot==1)','o','MarkerEdgeColor',colours(k,:),'MarkerFaceColor',colours(k,:),'MarkerSize',5); %plots data used in fitting
            plot(kap_plot,radii_plot(filter_plot==0)', phases_plot(filter_plot==0)','x','MarkerEdgeColor',colours(k,:),'MarkerSize',7); %plots excluded data
            errorbar(kap_plot,radii_plot,phases_plot,err_plot,'k','Marker','none', 'LineStyle','none', 'Color',colours(k,:))            
            
          
            plot(kap_plot,radii_plot',(phi(k,:))','-','Color',colours(k,:)); %plots curve fitted to individual data set
            % plot(radii_plot',phi_all(k,:),'--k','LineWidth',2); %plots curve fitted to all data for this period
            plot(kap_plot,calc_lag_r,calc_lag_phi,'--k','LineWidth',2); %plots curve fitted to all data for this period
        end

        %plot most significant data marker.
        plot(kap_plot,radii_mm(((most_significant(1)-1)*number_files )+ most_significant(2)), Phases_to_plot(((most_significant(1)-1)*number_files )+ most_significant(2))','o','MarkerEdgeColor','k','MarkerSize',8,'LineWidth',2); %plots data used in fitting
        
        %plots position of r0
        ylimits = [min(min(Phases_to_plot))-.05 max(max(Phases_to_plot))+.05];
        plot(kap_plot,[r0_mm r0_mm],ylimits,'k');

        set(kap_plot, 'NextPlot', 'replace'); %-- release figure without adding text to command window    
        
        %update display information
        set(info{1},'String',expt_name);
        set(info{2},'String',force(1));
        set(info{3},'String',temp(1));
        set(info{4},'String',which_set);
        set(radius_0_disp,'String', sprintf('%7.3f',r0_mm));
        set(kappa_display{1},'String', sprintf('%7.3f',kappa));
        for x = 1:number_files
            set(period_display{x+1},'String', sprintf('%7.3f',periods(x)));
            set(kappa_display{x+1},'String', sprintf('%7.3f',kappa_guesses(x)));
        end

        %wait for keyboard or mouse button to be pressed
        [x y button] = ginput(1);

        if isempty(x) == 1
            done = 1;
        elseif button == 115; %115 is 's'
            prompt = {'Offset between data sets'};
            dlg_title = 'Change offsets';
            num_lines = 1;
            def = {num2str(spacing)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if isempty(answer) == 0
                spacing = str2num(answer{1});
            end
        elseif button == 122 %122 is z
            zoom yon %allows zooming in the y direction on the figure

            waitfor(kap_plot,'Ylim') %allows one go at rescaling each time

            zoom off
            y_scaling = get(kap_plot,'YLim'); %defined to keep the rescaling subsequently
            
        elseif button == 120 || button == 88 % x or X.
            %             ending = questdlg...
            % here neeed a dlg asking if the data is rubbish and so end the
            % analysis of this data set. Then need to write to the output
            % file saying 'foil pair n at these conditions is no good!'
            
        elseif button == 102 || button == 70 % F or f - set r0 in solver.
            prompt = ['Value for r_0 (-1 to solve for r_0 otherwise #>0; maximum: ',num2str(max(max(radii_mm))),')'];
            dlg_title = 'R0 value';
            if R_solve == 1
                show = -1;
            else
                show = r0_mm;
            end
            def = {num2str(show)};
            answer = inputdlg(prompt,dlg_title,1,def);
            if isempty(answer) == 1
                %nothing in here - the if loop just catches the cancel from
                %the dialogue and stops the programme crashing
            elseif str2num(answer{1}) == -1
                R_solve = 1;
                r0_mm = 1;
            else
                R_solve = 0;
                r0_mm = str2num(answer{1});
            end
            
        elseif button >= 48 && button <= 48+length(kappa_guesses) %number between 0 and number data files
            change = button - 48;
            if change == 0 %replaces the global fit with the mean of the individualt fits.
                kappa = mean(kappa_guesses);
            else %replaces kappa for whichever period is chosen with that of the global fit.
                kappa_guesses(change) = kappa;
                phi0(change) = phi_nought(change);
            end

        elseif button == 116 || button == 84 %T or t -- sets up a template - which is saved later
            template = filter;
	    
        else
            %scale the distances so that using distances on screen. 
            x_scal = get(kap_plot,'Xlim');
            y_scal = get(kap_plot,'YLim');
            axes_pos = get(kap_plot,'Position');
            fig_pos = get(fig,'Position');
            axes_size = [axes_pos(3)*fig_pos(3) axes_pos(4)*fig_pos(4)];
            
            %position of chosen place in figure in pixels from corner of axes
            x_pix = (x - x_scal(1)) / (x_scal(2) - x_scal(1)) * axes_size(1);
            y_pix = (y - y_scal(1)) / (y_scal(2) - y_scal(1)) * axes_size(2);
            
            %position of radii and phases in figure in pixels from corner of axes
            rad_pix = (radii_mm - x_scal(1)) / (x_scal(2) - x_scal(1)) * axes_size(1);
            phases_pix = (Phases_to_plot - y_scal(1)) / (y_scal(2) - y_scal(1)) * axes_size(2);
            
            %distances between point and data in pixels
            x_dist = rad_pix - x_pix;
            y_dist = phases_pix - y_pix;
            
%             %distance between point and data in 'data-space'
%             x_dist = radii_mm - x;
%             y_dist = Phases_to_plot - y;

            distances =  sqrt(x_dist.^2 + y_dist.^2);
            [a b]= find(distances == min(min(distances)));

            if button == 1 % left mouse button - removes data point
                filter(a,b) = 0;
            elseif button == 3 % right mouse button - adds data point
                filter(a,b) = 1;
            elseif button == 112 || button == 88 % p or P  -- adds 2pi 
                phases(a,b) = phases(a,b) + 2*pi;
            elseif button == 109 || button == 77 %M or m  -- removes 2pi
                phases(a,b) = phases(a,b) - 2*pi;
            elseif button == 108 || button == 76 %L or l -- removes all data to left of mouse pointer.
                cull = x_dist < 0;
                filter(cull) = 0;
            elseif button == 114 || button == 82 %R or r -- removes all data to right of mouse pointer.
                cull = x_dist > 0;
                filter(cull) = 0;
            end
        end
    end

    %calculates average lengths between the boxes used to calcualte the
    %thermal diffusivity
    for d = 1 : number_files
        used = filter(d,:);
        length_to_use = lengths(d,:);
        mean_length(d) = mean(length_to_use(used == 1));
    end
    mean_length_all = mean(lengths(filter == 1));

    %% error calc
    try % the try statement allows for the catching of errors and allows the script to continue running...
        which_period = period(filter==1);
        good_omega = omegas(filter==1);
        phases_use = phases(filter==1);
        radii_use = radii_mm(filter==1);
        weights_use = weights(filter==1);
        
        % finds error on plus side of kappa by finding the minimum in acceleration in a misfit curve.
        % the half height between this residual and that at the best solution is
        % then used to find the errors.
        kappa_range_potential = (1:0.1:5) * kappa;
        
        kappa_range = kappa_range_potential(1:2);
        if R_solve == 1
            coeffs = [phi_nought r0_mm];
            values = kappa_range(1);
        else
            coeffs = phi_nought;
            values = [kappa_range(1) r0_mm];
        end

        %residual for kappa_range1
        residuals(1) = sqrt(misfit(coeffs,radii_use,phases_use,good_omega,which_period,weights_use,values));

        %residual for kapp_range2
        values(1) = kappa_range(2);
        best = fminsearch(@misfit,coeffs,[],radii_use,phases_use,good_omega,which_period,weights_use,values);
        residuals(2) = sqrt(misfit(best,radii_use,phases_use,good_omega,which_period,weights_use,values));

        slope_res(1) = (residuals(2) - residuals(1))/(kappa_range(2)-kappa_range(1));
        kappa_slope = (kappa_range(2)+kappa_range(1))/2;

        y = 3;
        x = 0;
        while x <= 2
            kappa_range = [kappa_range kappa_range_potential(y)];
            
            values(1) = kappa_range(y);
            best = fminsearch(@misfit,coeffs,[],radii_use,phases_use,good_omega,which_period,weights_use,values);
            residuals(y) = sqrt(misfit(best,radii_use,phases_use,good_omega,which_period,weights_use,values));

            slope_res(y-1) = (residuals(y) - residuals(y-1))/(kappa_range(y)-kappa_range(y-1));
            kappa_slope(y-1) = (kappa_range(y)+kappa_range(y-1))/2;

            curvature(y-2) = (slope_res(y-1) - slope_res(y-2))/(kappa_range(y-1)-kappa_range(y-2));
            kappa_accel(y-2) = (kappa_slope(y-1) + kappa_slope(y-2))/2;

            if y>=4
                if curvature(y-2) > curvature(y-3)
                    x = x+1;
                end
            end
            y = y + 1;
        end

        %difines kink in missfit curve
        max_curvature = find(min(curvature) == curvature);
        kappa_max_curve = kappa_range(max_curvature);

        mis_fit_at_max_curve = residuals(max_curvature);
        misfit_at_halfheight = (mis_fit_at_max_curve-min(residuals))/2 + min(residuals);


        diff_from_half_height = abs(residuals - misfit_at_halfheight);
        where = find(min(diff_from_half_height) == diff_from_half_height);

        smaller_gaps = kappa_range(where-1) : (kappa_range(where+1)-kappa_range(where-1))/20 : kappa_range(where+1);
        splined_residuals = spline(kappa_range(where-1:where+1),residuals(where-1:where+1),smaller_gaps);
        differences = abs(splined_residuals - misfit_at_halfheight);
        upper = smaller_gaps(find(min(differences)== differences));
        err_plus =  (upper - kappa)/sqrt(sum(sum(filter)));

        %finds error on minus side of kappa
        kappa_range_potential = (.995:-0.005:.005) * kappa;
        for x = 1:length(kappa_range_potential)
            values(1) = kappa_range_potential(x);
            best = fminsearch(@misfit,coeffs,[],radii_use,phases_use,good_omega,which_period,weights_use,values);
            residuals_small(x) = sqrt(misfit(best,radii_use,phases_use,good_omega,which_period,weights_use,values));
            if residuals_small(x) > misfit_at_halfheight
                break
            end
        end
        pos = length(residuals_small);
        lower = kappa_range_potential(pos) - (residuals_small(pos) - misfit_at_halfheight)/((residuals_small(pos) - residuals_small(pos-1))/(kappa_range_potential(pos) - kappa_range_potential(pos-1)));
        err_minus = (kappa - lower)/sqrt(sum(sum(filter)));
        
    catch
        warning('Kappa_solve:No_errors_calculated','There was an problem with the error calculation')
        
        err_minus = -1;
        err_plus = -1;
        
    end
%         figure(2)
%         clf
%         plot(kappa_range, residuals,'o-', kappa_range_potential(1:length(residuals_small)), residuals_small, 'o-', smaller_gaps,splined_residuals,'r-', [lower upper],[misfit_at_halfheight misfit_at_halfheight],'kx-')
%         xlabel('Kappa');
%         ylabel('Misfit');
%         legend('residuals', 'residuals small', 'splined residuals')
%     
%         figure(3)
%         plot(kappa_accel, curvature,'o-')
%     
%         pause
        
    R__nought = r0_mm;
    save(strcat(expt_name,'.mat'),'R_solve','R__nought', '-APPEND');
    if exist('template','var') == 1
        save(strcat(expt_name,'.mat'),'template', '-APPEND');
    end
    
    
    %% write outputfile
    where = cd;
    marker = find(where == filesep,1,'last');
    directory = where(marker+1:end);
    outfile_name = strcat(where,filesep,expt_name,'_thermal_diffusivity.txt');

    out_exist = exist(outfile_name,'file');
    outfile = fopen(outfile_name,'a');

    %makes preposition for outfile so that it is known how many times the
    %script has been run -- also used to separate the times when kappa solved
    %for singly and multpile times
    if out_exist == 0
        titles = {  '#';
            'File_name';
            'Load_(UStons)';
            'Temperature_(deg_C)';
            'nth_pair_from_middle';
            'Period_(s)';
            'r0_(mm)';
            'Kappa_(mm2s-1)_all_data';
            'err+_(mm2s-1)';
            'err-_(mm2s-1)';
            'mean_length_all_(mm)';
            'Phi0_(rad)_all_data';
%             'Reduced_ChiSq_all_data'; %<---new
            'Kappa_(mm2s-1)_separatly';
            'Phi0_(rad)_separately';
            'mean_length_per_period';
% %             'Reduced_ChiSq_separately'; %<---new
            'Filter';
            };
        fprintf(outfile,'%s',titles{1});
        for x = 2:length(titles)
            fprintf(outfile,',%s',titles{x});
        end
        fprintf(outfile,'\r');
        run_times = 1;
    else
        outfile_contents = importdata(outfile_name);
        number = length(outfile_contents);
        last_row = outfile_contents(number);
        last_row = char(last_row);
        separations = find(last_row ==',');
        initial = last_row(1:separations-1);
        run_times = str2num(initial) + 1;
    end

    for i = 1 : number_files
        fprintf(outfile,    '%3.0f, ',  run_times);
        fprintf(outfile,    '%s, ',     strtrim(file_names(i,:)));
        fprintf(outfile,    '%5.0f, ',  force(i));
        fprintf(outfile,    '%5.0f, ',  temp(i));
        
        fprintf(outfile,    '%s, ',     which_set);
        fprintf(outfile,    '%6.4f, ',  periods(i));
        fprintf(outfile,    '%6.4f, ',  r0_mm);
        
        fprintf(outfile,    '%7.4f, ',  kappa);
        fprintf(outfile,    '%7.4f, ',  err_plus);
        fprintf(outfile,    '%7.4f, ',  err_minus);
        fprintf(outfile,    '%3.3f, ',  mean_length_all);
        fprintf(outfile,    '%3.4f, ',  -phi_nought(i));
%         fprintf(outfile,    '%3.5f, ',  ChiSq_all);

        if sum(filter(i,:)) ~= 0
            fprintf(outfile,    '%7.4f, ',  kappa_guesses(i));
            fprintf(outfile,    '%3.4f, ',  -phi0(i));
            fprintf(outfile,    '%3.3f, ',  mean_length(i));
        else
            fprintf(outfile,    '%s, ',     'No_value');
            fprintf(outfile,    '%s, ',     'No_value');
            fprintf(outfile,    '%s, ',     'No_value');
        end
        
        fprintf(outfile,    '%s ',      ['[ ',num2str(filter(i,:)),' ]']);
        fprintf(outfile,'\r');
    end
    fclose('all');
end


end %kappa_solve


% #########################################################################
%% subfunctions

function out = kappaonly(coefs,X,Y,omega,limitsk,weights)
kappa = abs(coefs(1));
phi0 = coefs(2);

DIFF =  CalcPhase(X, omega, kappa) - (Y - phi0);
%CalcPhase calcualtes the Bessell function (in this case equaion 2 in Xu et
%al. (2004) and (Y- phi0) is the phase lags measured from the data. DIFF is
%the defference between the model and the data.

WTD_DIFF = weights.*DIFF.^2;
out2 = sum(sum(WTD_DIFF));
if kappa < limitsk(1)
    out = out2 * 10e10;
elseif kappa > limitsk(2)
    out = out2 * 10e10;
else
    out = out2;
end

end %kappa only


function out = KappaFind(coefs,X,Y,omega,which,limitsk,weights)
R_solve = coefs(1);
kappa = abs(coefs(2));
phi_nought = coefs(3:end)';

all_phi_0 = phi_nought(:,ones(1,size(which,2))) .* which;
all_phi_0 = all_phi_0(which == 1);

DIFF =  CalcPhase(abs(X), omega, kappa) - (Y - all_phi_0);
%CalcPhase calcualtes the Bessell function (in this case equaion 2 in Xu et
%al. (2004) and (Y- phi0) is the phase lags measured from the data. DIFF is
%the defference between the model and the data.

WTD_DIFF = weights.*DIFF.^2;
out1 = sum(sum(WTD_DIFF));
% if kappa < limitsk(1)
%     out2 = inf;
% elseif R_nought > limitsr(2)
%     out2 = inf;
% else
%     out2 = out1;
% end
if kappa > limitsk(2)
    out = out1 * 10e10;
else
    out = out1;
end

end %KappaFind

% function out = kappasearch(coefs,X,Y,omega,which,limitsr,limitsk)
% kappa = abs(coefs(1));
% R_nought = coefs(2);
% phi_nought = coefs(3:end)';
% 
% all_phi_0 = phi_nought(which);
% 
% DIFF =  CalcPhase(abs(X - R_nought), omega, kappa) - (Y - all_phi_0);
% %CalcPhase calcualtes the Bessell function (in this case equaion 2 in Xu et
% %al. (2004) and (Y- phi0) is the phase lags measured from the data. DIFF is
% %the defference between the model and the data.
% 
% SQ_DIFF = DIFF.^2;
% out1 = sum(sum(SQ_DIFF));
% if R_nought < limitsr(1)
%     out2 = out1 * 1E10;
% elseif R_nought > limitsr(2)
%     out2 = out1 * 1E10;
% else
%     out2 = out1;
% end
% if kappa > limitsk(2)
%     out = out2 * 10e10;
% else
%     out = out2;
% end
% 
% end %kappasearch

function out = kappasearchWEIGHTED(coefs,X,Y,omega,which,limitsr,limitsk,weights)
kappa = abs(coefs(1));
R_nought = coefs(2);
phi_nought = coefs(3:end)';

all_phi_0 = phi_nought(which);

DIFF =  CalcPhase(abs(X - R_nought), omega, kappa) - (Y - all_phi_0);
%CalcPhase calcualtes the Bessell function (in this case equaion 2 in Xu et
%al. (2004) and (Y- phi0) is the phase lags measured from the data. DIFF is
%the defference between the model and the data.

WTD_DIFF = weights.*DIFF.^2;
out1 = sum(sum(WTD_DIFF));
if R_nought < limitsr(1)
    out2 = out1 * 1E10;
elseif R_nought > limitsr(2)
    out2 = out1 * 1E10;
else
    out2 = out1;
end
if kappa > limitsk(2)
    out = out2 * 10e10;
else
    out = out2;
end
% fprintf(1,'%5.3f \r',out)
end %kappasearchWEIGHTED

function residual = misfit(coefs,X,Y,omega,which,weights,values)
kappa = values(1);
if length(values) == 2
    phi_nought = coefs';
    R_nought = values(2);
else
    phi_nought = coefs(1:end-1)';
    R_nought = coefs(length(coefs));
end
all_phi_0 = phi_nought(which);
DIFF =  CalcPhase(abs(X - R_nought), omega, kappa) - (Y - all_phi_0);

SQ_DIFF = weights.*DIFF.^2;
residual = sum(sum(SQ_DIFF));

end %residual

function Values = CalcPhase(R,Omega,Kappa)
%calculates the Bessel function for all instances they occur in the script
A = R .* ((Omega./Kappa).^(1/2));
[berA, beiA] = klvna(A);
Values = atan(beiA./berA);

end %phasecalc

function [ber,bei]=klvna(x)
%Edited from mklvna.m
% -- reduced so that it only calcualtes ber and bei rather than all the
% kelvin functions; the other functionals have been removed.

% mklvna.m is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).

%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)

eps=1e-15;
if(x == 0);
    ber=1;
    bei=0;
    return;
end;
x2 = 0.25 .* x .* x;
x4 = x2 .* x2;
% if(abs(x)< 10.0d0);
ber=ones(size(x)); %1;
r=ones(size(x)); %1;
for m = 1:60;
    r=-0.25.*r./(m.*m)./(2.*m-1).^2.*x4;
    ber=ber+r;
    if(abs(r)< abs(ber)*eps)
        break
    end
end;
bei=x2;
r=x2;
for m = 1:60;
    r=-0.25.*r./(m.*m)./(2.*m+1).^2.*x4;
    bei=bei+r;
    if (abs(r)< abs(bei)*eps)
        break
    end
end

end %klvna


function RedChi2 = ReducedChiSquared(obs,model,weights)
%copied from chi2test on mathworks website.

%Normalize x
N = length(obs);
x = obs - model;

% x = (x-mean(x))/std(x); %standardization
x = x/std(x);

xp = [-inf, -1.6:.4:1.6, inf]; %tested bins
xp = xp(:);
E = 0.5*N*diff(erfc(-xp/sqrt(2))); %expected frequency
S = histc(x, xp);

O = S(1:end-1); %%observed frequency
figure(2), plot(xp(2:end),E,'k-',xp(2:end),O,'k.'); figure(1)

Chi2=sum((E-O).^2./E); %statistics
d = (length(xp)-1) - 2; %%degrees of freedom
RedChi2 = Chi2/d; %%reduced Chi^2 value.

end %ReducedChiSquared

%% Versions:
% v 1.8 - 13 August 2014
%   - Changed reading of files so that the ReadSineFit is now used instead
%   of code within this file.
% v 1.7.2 - 30 August 2011
%   - Changed filter so that removes only data with largest errors as default
%   - Highlight most significant data point in fitting.
% v 1.7.1 - 10 June 2011
%   - pointed plotting at axes so that can add extra figures and not plot things in the wrong figure.
% v 1.7 - prior to 10 June 2011
%   - rewritten display generating code to separate it from the numberical analysis part of the code.
%   - cleaned up the code somewhat.
% v 1.5.0.3 (10th February 2011)
%   - fixed problem that generates new images each loop/time run and associated bugs
% v 1.5.0.2
%   - fixed error generated when no weights present in sin_fits file.
% 31-December-2010
%   - Started version numbering on the script. Started at v1.5 because if has not been totally
%   rewritten at any stage but has been substancailly added to
% Pre 31-December-2010:
%   - There is now a 'template' for the filter which can be passed between data sets
%   - Weighted fits for kappa using errors on phases.
%   - Converted all length units in to mm, no more m -- Kappa is now in mm²/s
%   - Checks to see that all the amplitudes are positive - if not corrects the phase (and amp)
%   - The distances used in removing data are now in pixels rather than distance in data-space.
% 21-June-2010:
%   - Phi0 is now the correct value and no longer pi out.
% 14-June-2010:
%   - This file incorporates most of the developments made to the script in
%   the past 9 months. Including:
%       + the ability to set R0 when solving for Kappa
% 30th Sept 2009.
%   - version 'old-1' 
%       + formal error routine
%       + 'experimental details' panel in window.
% - in archive as 'old-2' (NOT YET)
%       + made now script to read infor from file title and this script now
%               just calls that