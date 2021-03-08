function SSD_displacement_compare(file_list)

action = size(file_list);

close all

SSD_sine_ending = '_SSD_sine_fits.txt';
SSD_FITS_ending = '_SSD_FITS.mat';
single_sine_ending = '_sine_fits.txt';
position_ending = '_position_change.txt';

disp_FITS_ending = '_sine_fits_FITS.mat';
%open Zn02_27tons_25C_30s_001_SSD_FITS.mat


%% get experiment options
% - if they dont exist then it is all in vain
expt = FileTitleInformation(file_list{1});
varfilename = char(strcat(expt(1), '.mat'));

%list of variables need to get/generate.
filevariables = {'length_type', 'get_rid', 'symm_type', 'scaling', 'search'};

%reads the experiment analysis options file. -- the file is called 'expt_name'.mat
if exist(varfilename, 'file')
    vars = load(varfilename, filevariables{:});
end


%% read position change data files and other files
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
    data.(file_list{x}).length_change              = length_change_hold;
    data.(file_list{x}).image_time                 = image_time_hold;
    data.(file_list{x}).reference_length           = reference_length_hold;
    data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
    data.(file_list{x}).SSD                        = SSDdata;
    data.(file_list{x}).Disp                       = DispData;
        
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
number_pairs = size(length_change_hold, 2);

%number of foils in the data set.
number_fits = size(data.(file_list{x}).SSD.A_all,1);

%compare quality of fits
if 1
    
    for x = 1 : size(file_list)
        
        amp_err_rat(x,:) = data.(file_list{x}).Disp.output_values.amplitudes_err{1} ./data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        amp(x,:) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
        phase_err_rat(x,:) = data.(file_list{x}).Disp.output_values.phases_err{1} ./data.(file_list{x}).SSD.combined_fits.phases_err{1};
        ph(x,:) = data.(file_list{x}).SSD.combined_fits.phases_err{1};
        per(x,:) = data.(file_list{x}).SSD.combined_fits.period_err{1};
    end

    
    
    
    
    median(phase_err_rat(:))
    max(phase_err_rat(:))
    min(phase_err_rat(:))
    median(amp_err_rat(:))
    max(amp_err_rat(:))
    min(amp_err_rat(:))
    
    keyboard

end



% recreate sinusoidal variations from the fits. 
% labelled as strain in the structure soring all the data.
for x = 1 : size(file_list)
    
    times = data.(file_list{x}).image_time;
    
    for y = 1:number_pairs
        
        %calculate the sine dispacements for the non SSD data 
        
        %sine displacements from displament data
        amplitude = data.(file_list{x}).Disp.output_values.amplitudes{1}(y);
        period    = data.(file_list{x}).Disp.output_values.period{1};
        phase     = data.(file_list{x}).Disp.output_values.phases{1}(y);
        times     = data.(file_list{x}).Disp.image_time;
        
        %sine fit from displacements.
        data.(file_list{x}).Disp.Strain(y,:) = Phases_SineFunc(times, period, amplitude, phase);
        
%         StrainDisp(y,:)  = amplitude * sin(times./period(x,1)*2*pi + Phases(x,1));
        %StrainDisp(2,:) = amplitude(x,2) * sin(times./period(x,2)*2*pi + Phases(x,2));
              
        %model displacement from SSD data
        amplitude = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
        period    = data.(file_list{x}).SSD.combined_fits.period{1};
        phase     = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
        times     = data.(file_list{x}).SSD.image_time(y,:);
        
        %sine fit of the SSD data.
        data.(file_list{x}).SSD.Strain(y,:) = Phases_SineFunc(times, period, amplitude, phase);
                %strain1  = amplitude(x,1) * sin(times./period(x,1)*2*pi + Phases(x,1));
        %strain2 = amplitude(x,2) * sin(times./period(x,2)*2*pi + Phases(x,2));
        
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
        
        
        %as = polyvalnm(surf_coef, data.(file_list{x}).SSD.pos_array(:,:,z), repmat(All_times',search_dist*2+1,1))
        %df = polyvalnm(surf_diff, data.(file_list{x}).SSD.pos_array(:,:,z), repmat(All_times',search_dist*2+1,1));
        
        %polyvalnm_solve(surf_diff, 0, All_times')
        All_times   = data.(file_list{x}).image_time - data.(file_list{x}).image_time(1);
        data.(file_list{x}).SSD.surfmins(y,:) = polyvalnm_solve2(surf_diff, 0, data.(file_list{x}).SSD.image_time(y,:));
        
        %end
        
        %save the minimum lines if they do not exist in the data files.

        
    end
         %keyboard   
end

%% plot stuff

for x = 1 : size(file_list)
    
    if number_pairs <= 4
        
        
        %sort data for plotting.
        All_times   = data.(file_list{x}).image_time - data.(file_list{x}).image_time(1); 
        All_length  = data.(file_list{x}).length_change;
        
        %All_BackG   = data.(file_list{x}).Disp.BackGround; %in disp calculation bg is removed before sine is fit.
        
        Disp_times  = data.(file_list{x}).Disp.image_time;
        Disp_strain = data.(file_list{x}).Disp.Strain;
        Disp_backG  = data.(file_list{x}).Disp.BackGround';
        
        All_Fitted = ismember(All_times, Disp_times);
        
        SSD_times   = data.(file_list{x}).SSD.image_time;
        SSD_strain  = data.(file_list{x}).SSD.Strain;
        
        %plot fits vs. displacement.
        %displacements and fits
        
        %plot all foils and fits together.
        figure(1)
        for y = 1:number_pairs
            subplot(number_pairs,1,y), cla, hold on,
            
            %plot displacment data
            plot(All_times, All_length(:,y), ':b.')
            
            %plot disp fit
            plot(Disp_times, Disp_strain(y,:)+Disp_backG(y,:), '-r')
%             plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:), '-r')
            %plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:) + data.(file_list{x}).Disp.BackGround(:,y)', '-r')
            
            %plot SSD fit
            %plot(data.(file_list{x}).SSD.image_time(y,:), data.(file_list{x}).SSD.Strain(y,:), '-k')
            plot(SSD_times(y,:), SSD_strain(y,:)+data.(file_list{x}).SSD.surfmins(y,:)-data.(file_list{x}).SSD.surfmins(y+1,:), '-k')
            
            hold off
            if y == 1
                title(file_list{x}, 'Interpreter', 'none')
            end
        end
        
        
        %plot individual fits and residuals
        for y = 1:number_pairs
            
            figure
            subplot(2,1,1), hold on,
            
            %plot displacment data
            plot(All_times, All_length(:,y), 'b.', 'MarkerSize', 7)
            
            %plot disp fit
            plot(Disp_times, Disp_strain(y,:)+Disp_backG(y,:), '-r')
%             plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:), '-r')
            %plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:) + data.(file_list{x}).Disp.BackGround(:,y)', '-r')
            
            %plot SSD fit
            %plot(data.(file_list{x}).SSD.image_time(y,:), data.(file_list{x}).SSD.Strain(y,:), '-k')
            plot(SSD_times(y,:), SSD_strain(y,:)+data.(file_list{x}).SSD.surfmins(y,:)-data.(file_list{x}).SSD.surfmins(y+1,:), '-k')
            
            hold off
            
            ylabel('Length change (pixels)')
            xlabel('Time (s)')
            title('Fits')
            
            subplot(2,1,2), hold on,
            
            SSD_fit = SSD_strain(y,:)+data.(file_list{x}).SSD.surfmins(y,:)-data.(file_list{x}).SSD.surfmins(y+1,:);
            SSDresiduals = All_length(:,y) - SSD_fit';
            OldFitDiffs = Disp_strain(y,:)+Disp_backG(y,:) - SSD_fit;
            
            %plot displacment data
            plot(All_times, SSDresiduals, 'b.')
            
            %plot disp fit
            plot(Disp_times, OldFitDiffs, '-r')
%             plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:), '-r')
            %plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:) + data.(file_list{x}).Disp.BackGround(:,y)', '-r')
            
            %plot SSD fit
            %plot(data.(file_list{x}).SSD.image_time(y,:), data.(file_list{x}).SSD.Strain(y,:), '-k')
            plot(SSD_times(y,:), zeros(size(SSD_times(y,:))), '-k')
            
            ylabel('Residuals (pixels)')
            xlabel('Time (s)')
            title('Residuals')
            set(gca, 'Ylim', [-0.05, 0.05])
            
%             if y == 1
%                 title(file_list{x}, 'Interpreter', 'none')
%             end
        end
        
        
        
        
%         polyvalnm_diff(p,dim,n,varargin)
        
%         keyboard
        
%         polyvalnm_diff(polyvalnm_coef2mat(data.(file_list{x}).SSD.A_all(1,:), 6,3), 2,1)
%         polyvalnm_solve(polyvalnm_diff(polyvalnm_coef2mat(data.(file_list{x}).SSD.A_all(1,:), 6,3), 2,1), 0, All_times')
        
        if number_fits == 2
            
            e = 2;
            s = 1;
            %plot stress--strain
            figure(2), cla, hold on,
            plot(All_length(All_Fitted,e)-Disp_backG(e,:)', All_length(All_Fitted,s)-Disp_backG(s,:)', ':b.')
            %plot disp fit
            plot(Disp_strain(e,:), Disp_strain(s,:), '-r')
            %             plot(data.(file_list{x}).Disp.image_time, data.(file_list{x}).Disp.Strain(y,:) + data.(file_list{x}).Disp.BackGround(:,y)', '-r')
            %plot SSD fit
            plot(SSD_strain(e,:), SSD_strain(s,:), '-k')
            
        end
        
        
        
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