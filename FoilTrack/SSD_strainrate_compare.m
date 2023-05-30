function SSD_strainrate_compare(file_list, varargin)


global data

action = size(file_list);

close all

%set file endings
box_files = '_boxes';

SSD_FITS_ending = '_SSD_FITS.mat';
window_positionrate_ending = '_window_PositionChangeRate.txt';
window_lengthrate_ending = '_window_LengthChangeRate.txt';
SSDwindow_positionrate_ending = '_SSDwindow_PositionChangeRate.txt';
SSDwindow_lengthrate_ending = '_SSDwindow_LengthChangeRate.txt';

%SSD_sine_ending = '_SSD_sine_fits.txt';
%single_sine_ending = '_sine_fits.txt';
%position_ending = '_position_change.txt';
%SSD_position_ending = '_SSD_positions.txt';
%SSD_length_ending = '_SSDwindow_lengths.txt';
%disp_FITS_ending = '_sine_fits_FITS.mat';



detail_time_step = 0.2;

e = 2;
s = 1;

sample_names = [];

iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'samples'
            sample_names = varargin{iarg+1};
            iarg = iarg + 2;
        case 'step time'
            detail_time_step = varargin{iarg+1};
            iarg = iarg + 2;
        case 'anelastic'
            e = varargin{iarg+1};
            s = varargin{iarg+1};
            iarg = iarg + 3;
        
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end



if nargin >= 2
    sample_name = varargin{1};
else
    sample_name = 'sample';
end
if nargin >= 3
    standard_name = varargin{2};
else
    standard_name = 'standard';
end
if nargin >= 4
    standard_name = varargin{2};
else
    standard_name = 'sample';
    sample_name = 'standard';
end



%% get experiment options
% - if they dont exist then it is all in vain
if ischar(file_list)
    file_list = {file_list};
end
expt = FileTitleInformation(file_list{1});
varfilename = char(strcat(expt(1), '.mat'));

%list of variables need to get/generate.
filevariables = {'length_type', 'get_rid', 'symm_type', 'scaling', 'search'};

%reads the experiment analysis options file. -- the file is called 'expt_name'.mat
if exist(varfilename, 'file')
    vars = load(varfilename, filevariables{:});
end


%% read position change data files and other files
if isempty(data)
    disp('read the data')

    for x = 1:length(file_list)

        %get box reference positions
        if strcmpi(file_list{x}(end-2:end),"SSD")
            boxes = load([file_list{x}(1:end-4), '_boxes']);
        else
            boxes = load([file_list{x}, '_boxes']);
        end

        %import the position change-rates from raw displacement
        if exist([file_list{x}, window_positionrate_ending], 'file')
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, window_positionrate_ending]);

            data.(file_list{x}).disp_position.rate         = foil_change_data_x;
            data.(file_list{x}).disp_position.image_time   = time_stamps_x;
            %data.(file_list{x}).disp_length.reference_positions = box_positions_x;
            %data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
            data.(file_list{x}).number_images              = size(number_images_x,1);

            data.number_fits  = size(data.(file_list{1}).disp_position.rate, 2);
        end        
        %import the length change-rates from raw displacement
        if exist([file_list{x}, window_lengthrate_ending], 'file')
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, window_lengthrate_ending]);

            data.(file_list{x}).disp_length.rate           = foil_change_data_x;
            data.(file_list{x}).disp_length.image_time     = time_stamps_x;
            %data.(file_list{x}).disp_length.reference_positions = box_positions_x;
            %data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
            data.(file_list{x}).number_images              = size(number_images_x,1);

            data.number_pairs = size(data.(file_list{1}).disp_length.rate, 2)-1;
        end

        %import the SSD position change-rates from raw displacement
        if exist([file_list{x}, SSDwindow_positionrate_ending], 'file')
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, SSDwindow_positionrate_ending]);

            data.(file_list{x}).SSD_position.rate         = foil_change_data_x;
            data.(file_list{x}).SSD_position.image_time   = time_stamps_x;
            %data.(file_list{x}).disp_length.reference_positions = box_positions_x;
            %data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
            data.(file_list{x}).number_images              = size(number_images_x,1);
        end        
        %import the SSD length change-rates from raw displacement
        if exist([file_list{x}, SSDwindow_lengthrate_ending], 'file')
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, SSDwindow_lengthrate_ending]);

            data.(file_list{x}).SSD_length.rate           = foil_change_data_x;
            data.(file_list{x}).SSD_length.image_time     = time_stamps_x;
            %data.(file_list{x}).disp_length.reference_positions = box_positions_x;
            %data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
            data.(file_list{x}).number_images              = size(number_images_x,1);
        end

        if exist([file_list{x}, SSD_FITS_ending], 'file')
            %import SSD sine fits_FITS files.
            SSDdata = load([file_list{x}, SSD_FITS_ending]);
            data.(file_list{x}).SSD = SSDdata;
            data.(file_list{x}).used_images = size(SSDdata.image_time,2);
        else
            data.(file_list{x}).used_images = size(data.(file_list{x}).image_time);
        end

%         if exist([file_list{x}, SSD_position_ending], 'file')
%             %import SSD windowed positions
%             [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
%                 number_images_x] = ReadPositionChangeFile([file_list{x}, SSD_position_ending]);
%             data.(file_list{x}).SSDpositions.foilpositions = foil_change_data_x;
%             data.(file_list{x}).SSDpositions.time_stamps = time_stamps_x;
%             data.(file_list{x}).SSDpositions.box_positions = box_positions_x;
%             data.(file_list{x}).SSDpositions.number_boxes = number_boxes_x;
%             data.(file_list{x}).SSDpositions.number_images = number_images_x;
%         end
% 
%         if exist([file_list{x}, SSD_length_ending], 'file')
%             %import SSD windowed lengths
%             [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
%                 number_images_x] = ReadPositionChangeFile([file_list{x}, SSD_length_ending]);
%             data.(file_list{x}).SSDlengths.foilpositions = foil_change_data_x;
%             data.(file_list{x}).SSDlengths.time_stamps = time_stamps_x;
%             data.(file_list{x}).SSDlengths.box_positions = box_positions_x;
%             data.(file_list{x}).SSDlengths.number_boxes = number_boxes_x;
%             data.(file_list{x}).SSDlengths.number_images = number_images_x;
%         end
    end
end




%% recreate rate of length change from SSD fits (if present).
% labelled as strain in the structure storing all the data.

for x = 1 : size(file_list)
    disp(['Calculate displacement and length change-rates in ', file_list{x}])

    % windowed displcacement rates from disp data read from file 

    % windowed length change rates from disp data read from file 

    % windowed displcacement rates from SSD data read from file 

    % windowed length change rates from SSD data read from file 

    % displacment change from whole data set. 
    
    %calculate the sine dispacements and background for the SSD data

    %
    times2     = data.(file_list{x}).SSD.image_time(1,:);
    times_more2 = [times2(1):detail_time_step:times2(end)];
    data.(file_list{x}).SSD.times_more = times_more2;

    if isfield(data.(file_list{x}), 'SSD')
        surface_order = data.(file_list{x}).SSD.analysis_vars.polysurfcoeff;
        search_dist = data.(file_list{x}).SSD.analysis_vars.search;
        for y = 1:data.number_fits
            surf_coef = polyvalnm_coef2mat(data.(file_list{x}).SSD.A_all(y,:), surface_order(2), surface_order(1));
            surf_diff = polyvalnm_diff(surf_coef, 2,1);
            data.(file_list{x}).SSD.surf_diff(y,:,:) = surf_diff;
            surf_rate = polyvalnm_diff(surf_diff, 1, 1);
            data.(file_list{x}).SSD.surf_rate(y,:,:) = surf_rate;
    
            %polyvalnm(polyvalnm_coef2mat(A,n,m), pos_array(:,:,x), image_time))
            %[n m] = [2 6]
            
            %             coef =
            %
            %        NaN       NaN   -0.0000    0.0000   -0.0000    0.0001    0.0002
            %        NaN   -0.0000   -0.0000    0.0014    0.0044   -0.1616   -0.0363
            %     0.0028    0.0043   -0.3003   -0.3208   49.0174   19.9569    7.6299
            
            
            %polyvalnm_solve(surf_diff, 0, All_times')
            All_times        = data.(file_list{x}).SSD.image_time(1,:);
            All_times_more   = data.(file_list{x}).SSD.times_more - data.(file_list{x}).SSD.image_time(1,1);
            
            data.(file_list{x}).SSD.surfmins(y,:)      = polyvalnm_solve2(surf_diff, 0, All_times);%data.(file_list{x}).SSD.image_time(y,:));
            data.(file_list{x}).SSD.surfmins_more(y,:) = polyvalnm_solve2(surf_diff, 0, All_times_more);%data.(file_list{x}).SSD.image_time(y,:));
            if 0%isfield(data.(file_list{x}), 'SSD_length')
                %SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
                %data.(file_list{x}).SSD.surfmins_at_SSDlengths(y,:) = polyvalnm_solve2(surf_diff, 0, SSDlength_times);
                All_times_for_SSD_positions = data.(file_list{x}).SSDpositions.time_stamps - data.(file_list{x}).SSD.image_time(1,1);
                data.(file_list{x}).SSDpositions.surfminsrate_at_SSDlengths(y,:) = polyvalnm_solve2(surf_rate, 0, All_times_for_SSD_positions);
    
            end
            if 0%isfield(data.(file_list{x}), 'disp_length')
                %SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
                %data.(file_list{x}).SSD.surfmins_at_SSDlengths(y,:) = polyvalnm_solve2(surf_diff, 0, SSDlength_times);
                All_times_for_SSD_positions = data.(file_list{x}).SSDpositions.time_stamps - data.(file_list{x}).SSD.image_time(1,1);
                data.(file_list{x}).SSDpositions.surfminsrate_at_SSDlengths(y,:) = polyvalnm_solve2(surf_diff, 0, All_times_for_SSD_positions);
    
            end
            if 1

                times_all = data.(file_list{x}).disp_length.image_time - data.(file_list{x}).SSD.image_time0;
                data.(file_list{x}).SSDpositions.surfminsrate_at_SSDlengths(y,:) = polyvalnm_solve2(surf_rate, 0, times_all);
            end
    
            % check the differentiation is correct.
            if 0
    
                as = polyvalnm(surf_coef, data.(file_list{x}).SSD.pos_array(:,:,y), repmat(All_times,search_dist*2+1,1));
    
                xy_details = repmat((data.(file_list{x}).SSD.pos_array(1,1,y):.001:data.(file_list{x}).SSD.pos_array(end,1,y))',1,length(All_times));
                t_details = repmat(All_times,size(xy_details,1),1);
                as_detailed = polyvalnm(surf_coef, xy_details, t_details);
    
                [asdf, sdfg] = min(as_detailed);
                as_mins = xy_details(sdfg,1);
                df = polyvalnm(surf_diff, data.(file_list{x}).SSD.pos_array(:,:,y), repmat(All_times,search_dist*2+1,1));
    
                data.(file_list{x}).SSD.surfmins(y,:) = polyvalnm_solve2(surf_diff, 0, All_times);%data.(file_list{x}).SSD.image_time(y,:));
    
                plot(All_times, data.(file_list{x}).SSD.surfmins(y,:), 'b.-', All_times, as_mins, 'r.')
    
    
                %figure
                subplot(1,2,1)
                a = data.(file_list{x}).SSD.pos_array(:,:,y);
                b= repmat(All_times,search_dist*2+1,1);
                c= as;
                scatter3(a(:), b(:), c(:))
                subplot(1,2,2)
                a = data.(file_list{x}).SSD.pos_array(:,:,y),
                b= repmat(All_times,search_dist*2+1,1);
                c= df;
                scatter3(a(:), b(:), c(:))
                pause
            end
            %save the minimum lines if they do not exist in the data files.
    
        end

        %model displacement from SSD data
        for y = 1:data.number_pairs
            amplitude = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
            period    = data.(file_list{x}).SSD.combined_fits.period{1};
            phase     = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
            times2     = data.(file_list{x}).SSD.image_time(1,:);
            times_more2 = [round(times2(1),2):detail_time_step:round(times2(end),2)];
            if exist('times_more1','var') && length(times_more1) ~= length(times_more2)
                times_more2 = times_more1;
            end
            data.(file_list{x}).SSD.times_more = times_more2;

            %sine fit of the SSD data.
            data.(file_list{x}).SSD.sinusoidrate(y,:) = Phases_SineFunc(times2, period, amplitude/period*2*pi, phase+pi/2);
            data.(file_list{x}).SSD.sinusoidrate_more(y,:) = Phases_SineFunc(times_more2, period, amplitude/period*2*pi, phase+pi/2);
            data.(file_list{x}).SSD.phases(y,:) = rem(times2,period)/period*2*pi;

            if 1
                data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:) = Phases_SineFunc(data.(file_list{x}).SSD_length.image_time - data.(file_list{x}).SSD.image_time0, period, amplitude/period*2*pi, phase+pi/2);
                data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:) = Phases_SineFunc(data.(file_list{x}).disp_length.image_time - data.(file_list{x}).SSD.image_time0, period, amplitude/period*2*pi, phase+pi/2);
            end
            %strain1  = amplitude(x,1) * sin(times./period(x,1)*2*pi + Phases(x,1));
            %strain2 = amplitude(x,2) * sin(times./period(x,2)*2*pi + Phases(x,2));
        end
    end

% 
% 
%     %get reference lengths and changes in length
%     ref_positions = data.(file_list{x}).SSD.single_fits.reference_position{1};
%     change_position = data.(file_list{x}).SSD.surfmins(:,1) - data.(file_list{x}).SSD.surfmins(:,end);
%     
%     for y = 1:number_pairs
%         
%         period_all(x,y) = data.(file_list{x}).SSD.combined_fits.period{1};
%         period_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.period_err{1};
%         
%         phases_all(x,y) = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
%         phases_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.phases_err{1}(y);
%         
%         amplitudes_all(x,y) = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
%         amplitudes_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1}(y);
%         
%         ref_length(x,y) = ref_positions(y+1) - ref_positions(y);
% 
%         change_length(x,y) = change_position(y+1) - change_position(y);
%         if strcmpi(data.(file_list{x}).Disp.output_vars.length_type, 'pair') ==1
%             change_length(x,y) = change_position(number_pairs+y) - change_position(y);
%         elseif strcmpi(data.(file_list{x}).Disp.output_vars.length_type, 'anelastic') ==1
%             change_length(x,y) = change_position(y+1) - change_position(y);
%         else
%             error('not implemented')
%         end
% 
%         bulk_strain(x,y) = change_length(x,y)./ref_length(x,y);
%         elapsed_time(x) = times2(end)-times2(1);
%         bulk_strain_rate(x,y) = bulk_strain(x,y)./(times2(end)-times2(1));
%         
%         number_images(x) = data.(file_list{x}).used_images;
%     end
         %keyboard   
end


%% plot stuff
disp('plotting stuff')



%% plot strain-rates
for x = 1 : size(file_list)    
    if 1;%number_pairs <= 4
        
        
        %plot all foils and fits together.
        figure(1)
        clf
%         make_it_tight = true;
%         subplot2 = @(m,n,p) subtightplot(m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
        
        for y = 1:data.number_pairs
            ax(y) = subplot(data.number_pairs,1,y); cla, hold on,
            
            %plot displacment data
            plot(data.(file_list{x}).disp_length.image_time - data.(file_list{x}).SSD.image_time0, data.(file_list{x}).disp_length.rate(:,y), 'b.') %windowed rates
            hold on
            %plot(data.(file_list{x}).SSD_length.image_time - data.(file_list{x}).SSD.image_time0 , data.(file_list{x}).SSD_length.rate(:,y), 'r.') %windowed rates from SSD

            if isfield(data.(file_list{x}), 'SSD')
                plot(data.(file_list{x}).SSD.times_more-data.(file_list{x}).SSD.times_more(1), ...
                    data.(file_list{x}).SSD.sinusoidrate_more(y,:), 'k-')
            end


            if 0
                times_all = data.(file_list{x}).SSD_length.image_time - data.(file_list{x}).SSD.image_time0;
                %data.(file_list{x}).SSDpositions.surfminsrate_at_SSDlengths(y,:) = polyvalnm_solve2(surf_rate, 0, times_all);

                plot(times_all, ...
                    data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:), 'c.-')
                
            end
%             %plot disp fit
%             if exist('Disp_times_more', 'var')
%                 plot(Disp_times_more, Disp_lengths_more(y,:), '-r') %displacement fit
%             end
%             if exist('SSD_times_more', 'var')
%                 plot(SSD_times_more, SSD_length_more(y,:), '-k') %SSD minima fit
%             end
%             plot(SSD_times_more, SSD_backG_more(y,:), '--k') %SSD length background
            %plot(max(SSD_times_more)/2, SSD_ref_length(y), 'o') %SSD reference length
            ylabel('Length change (pixels/s)')
            xlabel('Time (s)')
            hold off
            limits(:,y) = get(ax(y), 'Ylim');
            
            box on
        end
        
        
        %set y-axes to be a nice ratio of each other
        if 1
            extent = limits(2,:) - limits(1,:);
            ratio = max(extent(:))./extent;
            rr = floor(ratio);
            rr_diffs = ratio-rr;
            for y = 1:data.number_pairs
                if rr_diffs(y) < 0
                    
                    %mean 
                    ax_mean = mean(limits(:,y));
                    %range
                    rng = max(extent(:));%./(rr(y)-1);
                    
                    ax_min = ax_mean-rng/2;
                    ax_max = ax_mean+rng/2;
                    
                    set(ax(y), 'Ylim', [ax_min, ax_max]);
                elseif rr_diffs(y) > 0
                    % expand the axes.
                    
                    %mean 
                    ax_mean = mean(limits(:,y));
                    %range
                    rng = max(extent(:));%./rr(y);
                    
                    ax_min = ax_mean-rng/2;
                    ax_max = ax_mean+rng/2;
                    
                    set(ax(y), 'Ylim', [ax_min, ax_max]);
                end
            end
        end

        %move figures to be close to each other. Remove space
        pos=get(ax,'position');
        bottom=pos{end}(2);
        top=pos{1}(2)+pos{1}(4);
        plotspace=top-bottom;
        for y = length(pos):-1:1
            pos{y}(2)=bottom+(length(pos)-y)*plotspace/length(pos);
            pos{y}(2)=bottom+(length(pos)-y)*pos{y}(4);
            %pos{y}(2) = bottom-pos{y}(4);
            set(ax(y),'position',pos{y});
        end
        for y = 1:length(pos)
            if rem(y,2) ~= 0
                set(ax(y),'YAxisLocation','right');
            end
            if y< length(pos)
               set(ax(y),'xticklabel',[]);
            end
        end
        
        %label the axes
        if 1
            for y = 1:length(ax)
                if length(sample_names) >= y
                    title_str = sample_names{y};
                else
                    title_str = ['sample ',num2str(y)];
                end
                tit(y)=title(ax(y),['  ',char(96+y),'. ',title_str]);
                tit(y).Position = [ ax(y).XLim(1), ax(y).YLim(2) ];
                tit(y).HorizontalAlignment = 'left';
                tit(y).VerticalAlignment = 'top';
            end
        end
            
        tightfig
        saveimage(['StrainRateFits_',file_list{x}])

        %plot individual fits and residuals
        % plot as tag group.
        if 0
            mainfig = figure(2);
            tabgroup = uitabgroup(mainfig, 'Position', [.01 .01 .98 .98]);
            
            for y = 1:data.number_pairs

                
                tab(y)=uitab(tabgroup,'Title', sprintf('Sample_%i', y));
                axes('parent',tab(y))
                                
                subplot(2,1,1), hold on,

                %plot displacment data
                plot(All_times, All_length(:,y), 'b.')
                %plot disp fits
                if exist('Disp_times_more', 'var')
                    plot(Disp_times_more, Disp_lengths_more(y,:), '-r') %displacement fit
                end
                if exist('SSD_times_more', 'var')
                    plot(SSD_times_more, SSD_length_more(y,:), '-k') %SSD minima fit
                end
                
                if isfield(data.(file_list{x}), 'SSDpositions')
                    SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
                    windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
                    plot(SSDlength_times', windowed_lengths, 'g.')
                end
                hold off
                ylabel('Length change (pixels)')
                xlabel('Time (s)')
                title('Fits')
                x_lims = get(gca,'Xlim');
                
                subplot(2,1,2), hold on,

                %plot residuals
                plot(All_times, All_length(:,y)-SSD_length(y,:)', 'b.')

                %plot disp fits
                if exist('Disp_times_more', 'var')
                    plot(Disp_times_more, Disp_lengths_more(y,:)-SSD_length_more(y,:), '-r')
                end
                if exist('SSD_times_more', 'var')
                    plot(SSD_times_more, SSD_length_more(y,:)-SSD_length_more(y,:), '-k')
                end
                
                if isfield(data.(file_list{x}), 'SSDpositions')
                    SSDlength_times  = data.(file_list{x}).SSDpositions.time_stamps;
                    windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
                    %all fit at times of windowed data
                    %background
                    SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);

                    plot(SSDlength_times', windowed_lengths'-SSD_all_ref, 'g.')

                end
                ylabel('Residuals (pixels)')
                xlabel('Time (s)')
                title('Residuals')
                set(gca, 'Xlim', x_lims);
            end
        end

        %% plot residuals for different fits (if applicable)
        figure
        for y = 1:data.number_pairs

            pSSD  = polyfit(data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).SSD_length.rate(:,y),1)
            pDisp = polyfit(data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).disp_length.rate(:,y),1)


            ax(y) = subplot(1,data.number_pairs,y); cla, hold on,
            plot(ax(y), data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).SSD_length.rate(:,y),  'r.')
            plot(ax(y), data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).disp_length.rate(:,y), 'b.')
            plot(ax(y),[min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))], ...
                [min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))], 'k-')

            plot(ax(y),[min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))], ...
                polyval(pSSD, [min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))]) , 'r--' )

            plot(ax(y),[min(data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:)),max(data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:))], ...
                polyval(pDisp, [min(data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:)),max(data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:))]) , 'b--' )

            ylabel('point-by-point rates (pixels/s)')
            xlabel('Rate from new algorithm (pixels/s)')
            title(ax(y),['  ',char(96+y),'. ',title_str]);

            disp('Correlation ')
            corrcoef(data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).SSD_length.rate(:,y))
            corrcoef(data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).disp_length.rate(:,y))


            %axis square
            %tightfig
        end
        saveimage(['RateResiduals_',file_list{x}])

        % plot differences
        figure
        for y = 1:data.number_pairs
            (data.(file_list{x}).SSD_length.rate(:,y)-data.(file_list{x}).disp_length.rate(:,y))./data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:)'
            ax(y) = subplot(1,data.number_pairs,y); cla, hold on,
            plot(ax(y), data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:), (data.(file_list{x}).SSD_length.rate(:,y)-data.(file_list{x}).disp_length.rate(:,y))./data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:)',  'r.')
            %plot(ax(y), data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:), data.(file_list{x}).disp_length.rate(:,y), 'b.')
            %plot(ax(y),[min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))], ...
            %    [min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))], 'k-')
            %plot(ax(y),[min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))], ...
            %    polyval(pSSD, [min(data.(file_list{x}).SSD.sinusoidrate_more(y,:)),max(data.(file_list{x}).SSD.sinusoidrate_more(y,:))]) , 'k--' )

            ylabel('delta rates (pixels/s)')
            xlabel('Rate from new algorithm (pixels/s)')
            title(ax(y),['  ',char(96+y),'. ',title_str]);
            %axis square
            %tightfig
        end
        saveimage(['RateResiduals2_',file_list{x}])



        figure
        for y = 1:data.number_pairs
            ax(y) = subplot(1,data.number_pairs,y); cla, hold on,

            resid1 = data.(file_list{x}).SSD_length.rate(:,y)' - data.(file_list{x}).SSD_length.sinusoidrate_at_SSDrates(y,:);
            resid2 = data.(file_list{x}).disp_length.rate(:,y)' - data.(file_list{x}).disp_length.sinusoidrate_at_SSDrates(y,:);
            disp(num2str(y))
            disp('SSD windowed residuals')
            disp(['mean ', num2str(mean(resid1))])
            disp(['std ', num2str(std(resid1))])
            disp('windowed point-by point residuals')
            disp(['mean ', num2str(mean(resid2))])
            disp(['std ', num2str(std(resid2))])
            plot(ax(y), sort(resid1), (1:length(resid1))/length(resid1)*100, 'r.')
            plot(ax(y), sort(resid2), (1:length(resid2))/length(resid2)*100, 'b.')
%             hold on
%             if isfield(data.(file_list{x}), 'SSDpositions')
%                 SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
%                 SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);
%                 %SSD_all_ref = SSD_all_ref+;
%                 windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
%                 resid_windowed = windowed_lengths'-SSD_all_ref;%SSD_all_ref-windowed_lengths'+mean(windowed_lengths);
% 
% 
%                 plot(ax(y), sort(resid_windowed), (1:length(resid_windowed))/length(resid_windowed)*100, 'g.')
% 
% 
%                 disp('windowed residuals')
%                 disp(['mean ', num2str(mean(resid_windowed))])
%                 disp(['std ', num2str(std(resid_windowed))])
% 
% 
%             end


            ylabel('Cumulative distribution (%)')
            xlabel('Residual')
            title('Residuals')
        end





        %% plot stress--strain
        if data.number_fits ==3
            
            
            figure(3); clf
            hold on,
            colormap('hsv')
            plot((All_length(All_Fitted,s)-SSD_backG(s,:)')/SSD_ref_length(s), (All_length(All_Fitted,e)-SSD_backG(e,:)')/SSD_ref_length(e), ':k')
            scatter((All_length(All_Fitted,s)-SSD_backG(s,:)')/SSD_ref_length(s), (All_length(All_Fitted,e)-SSD_backG(e,:)')/SSD_ref_length(e), 9, SSD_phases(s,:), 'filled');
            if exist('Disp_times_more', 'var')
                plot(Disp_sinusoid_more(s,:)/Disp_ref_length(s), Disp_sinusoid_more(e,:)/Disp_ref_length(e), '-r', 'LineWidth', 1)
            end
            if exist('SSD_times_more', 'var')
                plot(SSD_sinusoidrate_more(s,:)/SSD_ref_length(s), SSD_sinusoidrate_more(e,:)/SSD_ref_length(e), '-k', 'LineWidth', 1)
            end
            
            asdf = gca;
            set(gca, 'Clim', [0,2*pi]);
            ylabel(['Strain ', standard_name])
            xlabel(['Strain ',sample_name])
            c = colorbar('location','east outside');
            c.Label.String = 'Phase angle';
            
            
            %set axes to be a nice ratio of each other
            x_extent = get(gca, 'Xlim');
            y_extent = get(gca, 'Ylim');
            x_rng = x_extent(2) -x_extent(1);
            y_rng = y_extent(2) -y_extent(1);
            if x_rng < y_rng
                ratio = y_rng/x_rng;
                rr = floor(ratio);
                new_rng = y_rng/rr;
                ax_mean = mean(x_extent);
                ax_min = ax_mean-new_rng/2;
                ax_max = ax_mean+new_rng/2;
                    
                set(gca, 'Xlim', [ax_min, ax_max]);
                
            elseif x_rng > y_rng
                ratio = x_rng/y_rng;
                rr = floor(ratio);
                new_rng = x_rng/rr;
                ax_mean = mean(y_extent);
                ax_min = ax_mean-new_rng/2;
                ax_max = ax_mean+new_rng/2;
                    
                set(gca, 'Ylim', [ax_min, ax_max]);
            end
            %axis square
            tightfig
            saveimage(['Sin_sin_',file_list{x}])
        end
        
        %keyboard
        
        %         plot(times, strain1, '-r', time_stamps(1:number_images(x),x), length1(1:number_images(x),x), 'b.'),
        %
        %         title(file_list{x}, 'Interpreter', 'none')
        % %     subplot(2,3,3)
        % %     plot(reference_strain, length1(1:number_images(x),x), '.b'),
        %     subplot(2,3,4:6)
        %     plot(times, strain2, '-r', time_stamps(:,x), length2(:,x), 'b.')
        % %     subplot(2,3,6)
        % %     plot(sample_strain, length2(1:number_images(x),x), '-.b'),
        
        %pause
        
    else
        error('There are too many fits to plot here')
    end
    
end
    

return



function saveimage(name) 

    s = name;
    % could use savefigure here. It may be better but for now keep as is
    set(gcf,'PaperPositionMode','auto');
    %adjustpdfpage(gcf);
    %tightfig;
    
    %plname = strcat(s, '.jpg');
    %print(gcf, '-djpeg', plname);
    plname = strcat(s, '.pdf');
    print(gcf, '-dpdf', plname);
    %plname = strcat(s, '.eps');
    %print(gcf, '-depsc', plname);
    
    
