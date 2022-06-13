function SSD_displacement_compare(file_list, varargin)


global data

action = size(file_list);

close all

SSD_sine_ending = '_SSD_sine_fits.txt';
SSD_FITS_ending = '_SSD_FITS.mat';
single_sine_ending = '_sine_fits.txt';
position_ending = '_position_change.txt';

disp_FITS_ending = '_sine_fits_FITS.mat';
%open Zn02_27tons_25C_30s_001_SSD_FITS.mat

detail_time_step = 0.2;


e = 2;
s = 1;

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
    for x = 1:size(file_list,1)

        %import the postion change data
        [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
            number_images_x] = ReadPositionChangeFile([file_list{x}, position_ending]);

        %import displacement sine fits.
        DispData = load([file_list{x}, disp_FITS_ending]);

        %import SSD sine fits_FITS files.
        SSDdata = load([file_list{x}, SSD_FITS_ending]);

        %combine boxes. -- from raw displacement data.
        [length_change_hold, image_time_hold, reference_length_hold, horizontal_position_pixels_hold] = ...
        CombineBoxes(vars.length_type, 'displacements', foil_change_data_x, time_stamps_x, ...
                        box_positions_x(3:4,:)', box_positions_x(1:2,:)');

        %stash all the data in a structure.
        data.(file_list{x}).SSD                        = SSDdata;
        %data.(file_list{x}).Disp                       = DispData;

        data.(file_list{x}).length_change              = length_change_hold;
        data.(file_list{x}).image_time                 = image_time_hold;
        data.(file_list{x}).reference_length           = reference_length_hold;
        data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
        data.(file_list{x}).number_images              = size(number_images_x,1);
        
        data.(file_list{x}).used_images                = size(SSDdata.image_time,2);

    end
end

% %calculate length changes from position change data
% for x = 1:size(file_list,1)
%     foil_change_data(nSSD.umber_images(x):end,1:data_size(2),x) = NaN;
%     time_stamps(number_images(x):end,x) = NaN;
%     
%     length1(:,x) = foil_change_data(:,1,x) - foil_change_data(:,2,x);
%     length2(:,x) = foil_change_data(:,2,x) - foil_change_data(:,3,x);
% end



%time_stamps = time_stamps - repmat(time_stamps(1,:), size(length1,1),1);

% 
%FIX ME: We should check here that hte image times matach between the two data sets


% %% read SSD data files
% for x = 1:length(file_list)
%     SSD_files{x} = [file_list{x}, SSD_sine_ending];
% end
% qc = ReadSineFit('open', SSD_files{:});
% period = ReadSineFit('period');
% Period_err = ReadSineFit('period error');
% Phases  = ReadSineFit('phases');
% Phases_err = ReadSineFit('phases error');
% amplitude = ReadSineFit('amplitudes');
% amplitude_err = ReadSineFit('amplitudes error');
% ref_length  = ReadSineFit('lengths');
% qc = ReadSineFit('close');




%number of data sets from displacement length changes
number_pairs = size(data.(file_list{1}).length_change, 2);

%number of foils in the data set.
number_fits = size(data.(file_list{1}).SSD.A_all,1);

%compare quality of fits
if 1
    
    for x = 1 : size(file_list)
        
        amp_err_rat(x,:)   = data.(file_list{x}).Disp.output_values.amplitudes_err{1} ./data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        amp(x,:)           = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        phase_err_rat(x,:) = data.(file_list{x}).Disp.output_values.phases_err{1} ./data.(file_list{x}).SSD.combined_fits.phases_err{1};
        ph(x,:)            = data.(file_list{x}).SSD.combined_fits.phases_err{1};
        per(x,:)           = data.(file_list{x}).SSD.combined_fits.period_err{1};
        
        SSD_amp_list(x,:)       = data.(file_list{x}).SSD.combined_fits.amplitudes{1};
        SSD_amp_list_error(x,:) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        
    end

    median(phase_err_rat(:));
    max(phase_err_rat(:));
    min(phase_err_rat(:));
    median(amp_err_rat(:));
    max(amp_err_rat(:));
    min(amp_err_rat(:));
    
end

%compare quality of fits
if 1
    
    for x = 1 : size(file_list)
        
        amp_err_rat(x,:)   = data.(file_list{x}).Disp.output_values.amplitudes_err{1} ./data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        amps(x,:)           = data.(file_list{x}).SSD.single_fits.amplitudes_err{1};
        phase_err_rat(x,:) = data.(file_list{x}).Disp.output_values.phases_err{1} ./data.(file_list{x}).SSD.combined_fits.phases_err{1};
        phs(x,:)            = data.(file_list{x}).SSD.single_fits.phases_err{1};
        pers(x,:)           = data.(file_list{x}).SSD.single_fits.period_err{1};
        
        %SSD_amp_list(x,:)       = data.(file_list{x}).SSD.combined_fits.amplitudes{1};
        %SSD_amp_list_error(x,:) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        
    end
    if 0
        median(amps(:))
        max(amps(:))
        min(amps(:))
        median(phs(:)/pi*180)
        max(phs(:)/pi*180)
        min(phs(:)/pi*180)
        median(pers(:))
        max(pers(:))
        min(pers(:))
    end
end


% recreate sinusoidal variations from the fits. 
% labelled as strain in the structure soring all the data.
for x = 1 : size(file_list)
    
    %times = data.(file_list{x}).image_time;
    %times_more = [data.(file_list{x}).image_time(1):0.1:data.(file_list{x}).image_time(end)];
    
    for y = 1:number_pairs
        
        %calculate the sine dispacements for the non SSD data 
        if 0
            %sine displacements from displament data
            amplitude = data.(file_list{x}).Disp.output_values.amplitudes{1}(y);
            period    = data.(file_list{x}).Disp.output_values.period{1};
            phase     = data.(file_list{x}).Disp.output_values.phases{1}(y);
            times1     = data.(file_list{x}).Disp.image_time;
            times_more1 = [times1(1):detail_time_step:times1(end)];
            data.(file_list{x}).Disp.times_more = times_more1;

            %sine fit from displacements.
            data.(file_list{x}).Disp.sinusoid(y,:) = Phases_SineFunc(times1, period, amplitude, phase);
            data.(file_list{x}).Disp.sinusoid_more(y,:) = Phases_SineFunc(times_more1, period, amplitude, phase);
        end
        
%         StrainDisp(y,:)  = amplitude * sin(times./period(x,1)*2*pi + Phases(x,1));
        %StrainDisp(2,:) = amplitude(x,2) * sin(times./period(x,2)*2*pi + Phases(x,2));
        
        if 1
            %model displacement from SSD data
            amplitude = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
            period    = data.(file_list{x}).SSD.combined_fits.period{1};
            phase     = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
            times2     = data.(file_list{x}).SSD.image_time(y,:);
            times_more2 = [round(times2(1),2):detail_time_step:round(times2(end),2)];
            if exist('times_more1','var') && length(times_more1) ~= length(times_more2)
                times_more2 = times_more1;
            end
            data.(file_list{x}).SSD.times_more = times_more2;

            %sine fit of the SSD data.
            data.(file_list{x}).SSD.sinusoid(y,:) = Phases_SineFunc(times2, period, amplitude, phase);
            data.(file_list{x}).SSD.sinusoid_more(y,:) = Phases_SineFunc(times_more2, period, amplitude, phase);
            data.(file_list{x}).SSD.phases(y,:) = rem(times2,period)/period*2*pi;
                    %strain1  = amplitude(x,1) * sin(times./period(x,1)*2*pi + Phases(x,1));
            %strain2 = amplitude(x,2) * sin(times./period(x,2)*2*pi + Phases(x,2));
        end
    end
    
    %calculate the sine dispacements and background for the SSD data
    surface_order = data.(file_list{x}).SSD.analysis_vars.polysurfcoeff;
    search_dist = data.(file_list{x}).SSD.analysis_vars.search;
    
    for y = 1:number_fits
        
        %for z = 1:4
        surf_coef = polyvalnm_coef2mat(data.(file_list{x}).SSD.A_all(y,:), surface_order(2), surface_order(1));
        surf_diff = polyvalnm_diff(surf_coef, 2,1);
        
        %polyvalnm(polyvalnm_coef2mat(A,n,m), pos_array(:,:,x), image_time))
        %[n m] = [2 6]
        
        %             coef =
        %
        %        NaN       NaN   -0.0000    0.0000   -0.0000    0.0001    0.0002
        %        NaN   -0.0000   -0.0000    0.0014    0.0044   -0.1616   -0.0363
        %     0.0028    0.0043   -0.3003   -0.3208   49.0174   19.9569    7.6299
        
        
        %polyvalnm_solve(surf_diff, 0, All_times')
        All_times        = data.(file_list{x}).SSD.image_time(y,:) - data.(file_list{x}).SSD.image_time(y,1);
        All_times_more   = data.(file_list{x}).SSD.times_more - data.(file_list{x}).SSD.image_time(1,1);
        
        data.(file_list{x}).SSD.surfmins(y,:)      = polyvalnm_solve2(surf_diff, 0, All_times);%data.(file_list{x}).SSD.image_time(y,:));
        data.(file_list{x}).SSD.surfmins_more(y,:) = polyvalnm_solve2(surf_diff, 0, All_times_more);%data.(file_list{x}).SSD.image_time(y,:));
                
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
        end
                
        %save the minimum lines if they do not exist in the data files.

    end
    
    %get reference lengths and changes in length
    ref_positions = data.(file_list{x}).SSD.single_fits.reference_position{1};
    change_position = data.(file_list{x}).SSD.surfmins(:,1) - data.(file_list{x}).SSD.surfmins(:,end);
    
    for y = 1:number_pairs
        
        period_all(x,y) = data.(file_list{x}).SSD.combined_fits.period{1};
        period_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.period_err{1};
        
        phases_all(x,y) = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
        phases_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.phases_err{1}(y);
        
        amplitudes_all(x,y) = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
        amplitudes_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1}(y);
        
        ref_length(x,y) = ref_positions(y+1) - ref_positions(y);
        change_length(x,y) = change_position(y+1) - change_position(y);
        
        bulk_strain(x,y) = change_length(x,y)./ref_length(x,y);
        elapsed_time(x) = times2(end)-times2(1);
        bulk_strain_rate(x,y) = bulk_strain(x,y)./(times2(end)-times2(1));
        
        number_images(x) = data.(file_list{x}).used_images;
    end
         %keyboard   
end



keyboard
%% plot stuff

for x = 1 : size(file_list)
    
    if number_pairs <= 4
        
        
        %sort data for plotting.
        All_times   = data.(file_list{x}).image_time - data.(file_list{x}).image_time(1); 
        All_length  = data.(file_list{x}).length_change + data.(file_list{x}).reference_length;
        
        %All_BackG   = data.(file_list{x}).Disp.BackGround; %in disp calculation bg is removed before sine is fit.
        
        Disp_times         = data.(file_list{x}).Disp.image_time;
        Disp_times_more    = data.(file_list{x}).Disp.times_more;
        Disp_sinusoid      = data.(file_list{x}).Disp.sinusoid;
        Disp_sinusoid_more = data.(file_list{x}).Disp.sinusoid_more;
        
        Disp_ref_length    = data.(file_list{x}).Disp.output_values.reference_length{1};
        
        Disp_backG         = data.(file_list{x}).Disp.BackGround'+Disp_ref_length';
        Disp_backG_more = [];
        for y = 1:size(Disp_backG,1)
            Disp_backG_more(y,:) = interp1(Disp_times, Disp_backG(y,:), data.(file_list{x}).Disp.times_more);
        end
        
        Disp_lengths       = Disp_sinusoid +Disp_backG;
        Disp_lengths_more  = Disp_sinusoid_more+Disp_backG_more;
        
        All_Fitted = ismember(All_times, Disp_times);
        
        SSD_times              = data.(file_list{x}).SSD.image_time;
        SSD_times_more         = data.(file_list{x}).SSD.times_more;%image_time;
        SSD_sinusoid      = data.(file_list{x}).SSD.sinusoid;
        SSD_sinusoid_more = data.(file_list{x}).SSD.sinusoid_more;%sinusoid;
        
        SSD_ref_length    = data.(file_list{x}).SSD.combined_fits.reference_length{1};
        
        SSD_backG         = data.(file_list{x}).SSD.surfmins(1:end-1,:)-data.(file_list{x}).SSD.surfmins(2:end,:)+ round(data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end))-round(data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1));
        SSD_backG_more    = data.(file_list{x}).SSD.surfmins_more(1:end-1,:)-data.(file_list{x}).SSD.surfmins_more(2:end,:)+ round(data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end))-round(data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1));
                
        SSD_length        = SSD_sinusoid + SSD_backG;
        SSD_length_more   = SSD_backG_more + SSD_sinusoid_more;
        
        SSD_phases        = mod(data.(file_list{x}).SSD.phases + data.(file_list{x}).SSD.combined_fits.phases{1}', 2*pi);

        
        
        %plot all foils and fits together.
        figure(1)
        clf
%         make_it_tight = true;
%         subplot2 = @(m,n,p) subtightplot(m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
        
        for y = 1:number_pairs
            ax(y) = subplot(number_pairs,1,y); cla, hold on,
            
            %plot displacment data
            plot(All_times, All_length(:,y), 'b.') %SSD minimia
            
            %plot disp fit
            plot(Disp_times_more, Disp_lengths_more(y,:), '-r') %displacement fit
            plot(SSD_times_more, SSD_length_more(y,:), '-k') %SSD minima fit
            plot(SSD_times_more, SSD_backG_more(y,:), '--k') %SSD length background
            %plot(max(SSD_times_more)/2, SSD_ref_length(y), 'o') %SSD reference length
            ylabel('Length (pixels)')
            xlabel('Time (s)')
            hold off
            limits(:,y) = get(ax(y), 'Ylim');
            
            box on
            
        end
        
        %set y-axes to be a nice ratio of each other
        extent = limits(2,:) - limits(1,:);
        ratio = max(extent(:))./extent;
        rr = floor(ratio);
        rr_diffs = ratio-rr;
        for y = 1:number_pairs
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
        
        %move figures to be close to each other. Remove space
        
        pos=get(ax,'position');
        for y = 2:length(pos)
            bottom=pos{y-1}(2);
            top=pos{y}(2)+pos{y}(4);
            plotspace=top-bottom;
            pos{y}(2)=bottom+plotspace/2;
            pos{y}(2) = bottom-pos{y}(4);
            set(ax(y),'position',pos{2});
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
                if y==s
                    title_str = sample_name;
                elseif y==e
                    title_str = standard_name;
                else
                    title_str = num2str(y);
                end
                tit(y)=title(ax(y),['  ',char(96+y),'. ',title_str]);
                tit(y).Position = [ ax(y).XLim(1), ax(y).YLim(2) ];
                tit(y).HorizontalAlignment = 'left';
                tit(y).VerticalAlignment = 'top';
            end
        end
            
        
        
        
        tightfig
        saveimage(['Fits_',file_list{x}])
        %keyboard
        
        %plot individual fits and residuals
        % plot as tag group.
        if 1
            mainfig = figure(2);
            tabgroup = uitabgroup(mainfig, 'Position', [.01 .01 .98 .98]);
            
            for y = 1:number_pairs

                
                tab(y)=uitab(tabgroup,'Title', sprintf('Sample_%i', y));
                axes('parent',tab(y))
                                
                subplot(2,1,1), hold on,

                %plot displacment data
                plot(All_times, All_length(:,y), 'b.')
                %plot disp fits
                plot(Disp_times_more, Disp_lengths_more(y,:), '-r')
                plot(SSD_times_more, SSD_length_more(y,:), '-k')
                hold off
                ylabel('Length change (pixels)')
                xlabel('Time (s)')
                title('Fits')
                x_lims = get(gca,'Xlim');
                
                subplot(2,1,2), hold on,

                %plot displacment data
                plot(All_times(ismember(round(All_times,2),round(SSD_times(1,:),2))), All_length(ismember(round(All_times,2),round(SSD_times(1,:),2)),y)-SSD_length(y,:)', 'b.')

                %plot disp fits
                plot(Disp_times_more, Disp_lengths_more(y,:)-SSD_length_more(y,:), '-r')
                plot(SSD_times_more, SSD_length_more(y,:)-SSD_length_more(y,:), '-k')
                
                ylabel('Residuals (pixels)')
                xlabel('Time (s)')
                title('Residuals')
                set(gca, 'Xlim', x_lims);
                %set(gca, 'Ylim', [-0.05, 0.05])

    %             if y == 1
    %                 title(file_list{x}, 'Interpreter', 'none')
    %             end
            end
        end
        
        if number_fits == 3
            
            %plot stress--strain
            figure(3); clf
                        
            plot((All_length(All_Fitted,s)-SSD_backG(s,:)')/SSD_ref_length(s), (All_length(All_Fitted,e)-SSD_backG(e,:)')/SSD_ref_length(e), ':k')
            hold on,
            colormap('hsv')
            scatter((All_length(All_Fitted,s)-SSD_backG(s,:)')/SSD_ref_length(s), (All_length(All_Fitted,e)-SSD_backG(e,:)')/SSD_ref_length(e), 9, SSD_phases(s,:), 'filled');
            
            plot(Disp_sinusoid_more(s,:)/Disp_ref_length(s), Disp_sinusoid_more(e,:)/Disp_ref_length(e), '-r', 'LineWidth', 1)
            plot(SSD_sinusoid_more(s,:)/SSD_ref_length(s), SSD_sinusoid_more(e,:)/SSD_ref_length(e), '-k', 'LineWidth', 1)
            
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

    for x = 1:size(file_list,1)
    
    [maxwell_E(x), maxwell_nu(x)] = compute_maxwell_model(Phases(x,1), amplitude(x,1), ...
    Phases(x,2), amplitude(x,2), 350, period(x,1))

end
plot(real(maxwell_E))
figure
plot(maxwell_nu)




function saveimage(name) 

    s = name;
    % could use savefigure here. It may be better but for now keep as is
    set(gcf,'PaperPositionMode','auto');
    adjustpdfpage(gcf);
    %tightfig;
    
    %plname = strcat(s, '.jpg');
    %print(gcf, '-djpeg', plname);
    plname = strcat(s, '.pdf');
    print(gcf, '-dpdf', plname);
    %plname = strcat(s, '.eps');
    %print(gcf, '-depsc', plname);
    
    
