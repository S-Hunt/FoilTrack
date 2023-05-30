function SSD_displacement_compare(file_list, varargin)


global data

action = size(file_list);

close all

%set file endings
box_files_ending = '_boxes';

disp_position_ending = '_position_change.txt';
disp_sine_single_ending = '_sine_fits_single.txt';
disp_sine_ending = '_sine_fits.txt';
disp_FITS_ending = '_sine_fits_FITS.mat';

SSD_sine_ending = '_SSD_sine_fits.txt';
SSD_FITS_ending = '_SSD_FITS.mat';
SSD_position_ending = '_SSD_positions.txt';
SSD_length_ending = '_SSD_lengths.txt';

window=7;
%open Zn02_27tons_25C_30s_001_SSD_FITS.mat

detail_time_step = 0.2;


%for anelastic experiments
%set the reference and the sample
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



%% make sure the file list is good.
for x = 1:length(file_list)
    if strcmpi(file_list{x}(end-2:end),"SSD")
        file_list{x} = file_list{x}(end-3:end);
    end
end

%% read position fit files and other files
if isempty(data)
    disp('read the data')

    for x = 1:length(file_list)

        %get box reference positions
        boxes = load([file_list{x}, '_boxes']);

        %import the displacement (non SSD) data
        if exist([file_list{x}, disp_position_ending], 'file') == 2
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, disp_position_ending]);

            %fix reference lengths from rounded values in file
            box_positions_x=[boxes.boxY';boxes.boxX'];
                        
            %combine boxes. -- from raw displacement data.
            [length_change_hold, image_time_hold, reference_length_hold, horizontal_position_pixels_hold] = ...
            CombineBoxes(vars.length_type, 'displacements', foil_change_data_x, time_stamps_x, ...
                            box_positions_x(3:4,:)', box_positions_x(1:2,:)');
            
            %store data in a structure.
            data.(file_list{x}).image_time                 = time_stamps_x;
            data.(file_list{x}).reference_positions        = box_positions_x;
            data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
            data.(file_list{x}).number_images              = size(number_images_x,1);

            data.(file_list{x}).Disp.position_change            = foil_change_data_x;
            data.(file_list{x}).Disp.position = data.(file_list{x}).Disp.position_change + mean(data.(file_list{x}).reference_positions(1:2,:));
            data.(file_list{x}).Disp.length_change              = length_change_hold;
            data.(file_list{x}).Disp.reference_length           = reference_length_hold;
        end
% 
%         %import the displacement (non SSD) data
%         if exist([file_list{x}, disp_position_ending], 'file') == 2
%             [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
%                 number_images_x] = ReadPositionChangeFile([file_list{x}, disp_position_ending]);
% 
%             %fix reference lengths from rounded values in file
%             box_positions_x=[boxes.boxY';boxes.boxX'];
%                         
%             %combine boxes. -- from raw displacement data.
%             [length_change_hold, image_time_hold, reference_length_hold, horizontal_position_pixels_hold] = ...
%             CombineBoxes(vars.length_type, 'displacements', foil_change_data_x, time_stamps_x, ...
%                             box_positions_x(3:4,:)', box_positions_x(1:2,:)');
%             
%             %store data in a structure.
%             data.(file_list{x}).image_time                 = time_stamps_x;
%             data.(file_list{x}).reference_positions        = box_positions_x;
%             data.(file_list{x}).horizontal_position_pixels = horizontal_position_pixels_hold;
%             data.(file_list{x}).number_images              = size(number_images_x,1);
% 
%             data.(file_list{x}).Disp.position_change            = foil_change_data_x;
%             data.(file_list{x}).Disp.position = data.(file_list{x}).Disp.position_change + mean(data.(file_list{x}).reference_positions(1:2,:));
%             data.(file_list{x}).Disp.length_change              = length_change_hold;
%             data.(file_list{x}).Disp.reference_length           = reference_length_hold;
%         end



        %if all the displacement exist as a mat file them import them.
        if exist([file_list{x}, disp_sine_single_ending], 'file')==2
            %import single sine fits.
            ReadSineFit('open',[file_list{x}, disp_sine_single_ending]);
            data.(file_list{x}).Disp.SingleSineFits.period = ReadSineFit('Period');
            data.(file_list{x}).Disp.SingleSineFits.phases = ReadSineFit('Phases');
            data.(file_list{x}).Disp.SingleSineFits.phases_err = ReadSineFit('Phases error');
            data.(file_list{x}).Disp.SingleSineFits.amplitudes = ReadSineFit('Amplitudes');
            data.(file_list{x}).Disp.SingleSineFits.amplitudes_err = ReadSineFit('Amplitudes error');

            data.(file_list{x}).Disp.SingleSineFits.backgrounds = ReadSineFit('polynomials')
        end

        % import displacement sine fits from mat file
        if exist([file_list{x}, disp_FITS_ending], 'file')
            DispData = load([file_list{x}, disp_FITS_ending]);
            data.(file_list{x}).DispData = DispData;
            if exist('boxes','var')
                %fix reference lengths from rounded values in file
                data.(file_list{x}).Disp.output_values.reference_length{1} = reference_length_hold;
            end
        end
        %import SSD sine fits_FITS files.
        if exist([file_list{x}, SSD_FITS_ending], 'file')
            SSDdata = load([file_list{x}, SSD_FITS_ending]);
            data.(file_list{x}).SSD = SSDdata;
            data.(file_list{x}).used_images = size(SSDdata.image_time,2);
        else
            data.(file_list{x}).used_images = size(data.(file_list{x}).image_time);
        end

        %import SSD windowed positions
        if exist([file_list{x}, SSD_position_ending], 'file')
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, SSD_position_ending]);
            data.(file_list{x}).SSDpositions.foilpositions = foil_change_data_x;
            data.(file_list{x}).SSDpositions.time_stamps = time_stamps_x;
            data.(file_list{x}).SSDpositions.box_positions = box_positions_x;
            data.(file_list{x}).SSDpositions.number_boxes = number_boxes_x;
            data.(file_list{x}).SSDpositions.number_images = number_images_x;
        end

        %import SSD windowed lengths
        if exist([file_list{x}, SSD_length_ending], 'file')
            [foil_change_data_x, time_stamps_x, box_positions_x, number_boxes_x, ...
                number_images_x] = ReadPositionChangeFile([file_list{x}, SSD_length_ending]);
            data.(file_list{x}).SSDlengths.foilpositions = foil_change_data_x;
            data.(file_list{x}).SSDlengths.time_stamps = time_stamps_x;
            data.(file_list{x}).SSDlengths.box_positions = box_positions_x;
            data.(file_list{x}).SSDlengths.number_boxes = number_boxes_x;
            data.(file_list{x}).SSDlengths.number_images = number_images_x;
        end
    end
end

%number of data sets from displacement length changes
number_fits  = size(data.(file_list{1}).Disp.position_change, 2);
number_pairs = size(data.(file_list{1}).Disp.length_change, 2);

%number of foils in the data set.
%number_fits = size(data.(file_list{1}).SSD.A_all,1);



%% compare quality of fits
% for x = 1 : size(file_list)
%     if isfield(data.(file_list{x}), 'Disp') & isfield(data.(file_list{x}), 'SSD')
%         amp_err_rat(x,:)    = data.(file_list{x}).DispData.output_values.amplitudes_err{1} ./data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
%         amps(x,:)           = data.(file_list{x}).SSD.single_fits.amplitudes_err{1};
%         phase_err_rat(x,:)  = data.(file_list{x}).DispData.output_values.phases_err{1} ./data.(file_list{x}).SSD.combined_fits.phases_err{1};
%         phs(x,:)            = data.(file_list{x}).SSD.single_fits.phases_err{1};
%         pers(x,:)           = data.(file_list{x}).SSD.single_fits.period_err{1};
%         
%         %SSD_amp_list(x,:)       = data.(file_list{x}).SSD.combined_fits.amplitudes{1};
%         %SSD_amp_list_error(x,:) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1};
%     else
%         amp_err_rat(x,:)    = NaN;
%         amps(x,:)           = NaN;
%         phase_err_rat(x,:)  = NaN;
%         phs(x,:)            = NaN;
%         pers(x,:)           = NaN;
%     end
% end
% if 0
%     median(amps(:))
%     max(amps(:))
%     min(amps(:))
%     median(phs(:)/pi*180)
%     max(phs(:)/pi*180)
%     min(phs(:)/pi*180)
%     median(pers(:))
%     max(pers(:))
%     min(pers(:))
% end


%% recreate length changes and sinusoidal variations from the fits.
% labelled as strain in the structure storing all the data.

for x = 1 : size(file_list)
    disp(['Calculate displacements and length changes in ', file_list{x}])

    if isfield(data.(file_list{x}), 'Disp')
        %calculate the sine dispacements for the non SSD data
        for y = 1:number_pairs
            %sine displacements from displament data
            amplitude = data.(file_list{x}).DispData.output_values.amplitudes{1}(y);
            period    = data.(file_list{x}).DispData.output_values.period{1};
            phase     = data.(file_list{x}).DispData.output_values.phases{1}(y);
            times1     = data.(file_list{x}).DispData.image_time;
            times_more1 = [times1(1):detail_time_step:times1(end)];
            data.(file_list{x}).Disp.times_more = times_more1;

            %sine fit from displacements.
            data.(file_list{x}).Disp.sinusoid(y,:) = Phases_SineFunc(times1, period, amplitude, phase);
            data.(file_list{x}).Disp.sinusoid_more(y,:) = Phases_SineFunc(times_more1, period, amplitude, phase);
        end

        
        for y = 1:number_pairs
            %sine displacements from displament data
            amplitude = data.(file_list{x}).DispData.output_values.amplitudes{1}(y);
            period    = data.(file_list{x}).DispData.output_values.period{1};
            phase     = data.(file_list{x}).DispData.output_values.phases{1}(y);
            times1     = data.(file_list{x}).DispData.image_time;
            times_more1 = [times1(1):detail_time_step:times1(end)];
            data.(file_list{x}).Disp.times_more = times_more1;

            %sine fit from displacements.
            data.(file_list{x}).Disp.sinusoid(y,:) = Phases_SineFunc(times1, period, amplitude, phase);
            data.(file_list{x}).Disp.sinusoid_more(y,:) = Phases_SineFunc(times_more1, period, amplitude, phase);
        end
    end

    if isfield(data.(file_list{x}), 'SSD')

        %calculate the sine dispacements and background for the SSD data
        surface_order = data.(file_list{x}).SSD.analysis_vars.polysurfcoeff;
        search_dist = data.(file_list{x}).SSD.analysis_vars.search;
        for y = 1:number_fits

            %calcuate sinusoid length changes for each area of interest
            amplitude = data.(file_list{x}).SSD.single_fits.amplitudes{1}(y);
            period    = data.(file_list{x}).SSD.single_fits.period{1};
            phase     = data.(file_list{x}).SSD.single_fits.phases{1}(y);
            times2     = data.(file_list{x}).SSD.image_time(1,:);
            times_more2 = [round(times2(1),2):detail_time_step:round(times2(end),2)];
            if exist('times_more1','var') && length(times_more1) ~= length(times_more2)
                times_more2 = times_more1;
            end
            data.(file_list{x}).SSD.times_more = times_more2;

            data.(file_list{x}).SSD.single_sinusoid(y,:) = Phases_SineFunc(times2, period, amplitude, phase);
            data.(file_list{x}).SSD.single_sinusoid_more(y,:) = Phases_SineFunc(times_more2, period, amplitude, phase);
            data.(file_list{x}).SSD.single_phases(y,:) = rem(times2,period)/period*2*pi;


            %reproduce surface minima for SSD fits
    
            surf_coef = polyvalnm_coef2mat(data.(file_list{x}).SSD.A_all(y,:), surface_order(2), surface_order(1));
            surf_diff = polyvalnm_diff(surf_coef, 2,1);
            data.(file_list{x}).SSD.surf_diff(y,:,:) = surf_diff;
            %polyvalnm(polyvalnm_coef2mat(A,n,m), pos_array(:,:,x), image_time))
            %[n m] = [2 6]
            
            %             coef =
            %
            %        NaN       NaN   -0.0000    0.0000   -0.0000    0.0001    0.0002
            %        NaN   -0.0000   -0.0000    0.0014    0.0044   -0.1616   -0.0363
            %     0.0028    0.0043   -0.3003   -0.3208   49.0174   19.9569    7.6299
            

            All_times        = data.(file_list{x}).SSD.image_time(1,:) - data.(file_list{x}).SSD.image_time(1,1);
            All_times_more   = data.(file_list{x}).SSD.times_more - data.(file_list{x}).SSD.image_time(1,1);
            data.(file_list{x}).SSD.single_surfmins(y,:)      = polyvalnm_solve2(surf_diff, 0, All_times);%data.(file_list{x}).SSD.image_time(y,:));
            data.(file_list{x}).SSD.single_surfmins_more(y,:) = polyvalnm_solve2(surf_diff, 0, All_times_more);%data.(file_list{x}).SSD.image_time(y,:));

    
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

        %model sin component of length changes from SSD data
        for y = 1:number_pairs
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
            data.(file_list{x}).SSD.combined_sinusoid(y,:) = Phases_SineFunc(times2, period, amplitude, phase);
            data.(file_list{x}).SSD.combined_sinusoid_more(y,:) = Phases_SineFunc(times_more2, period, amplitude, phase);
            data.(file_list{x}).SSD.combined_phases(y,:) = rem(times2,period)/period*2*pi;
        end
        

        % calculate compelete model at times of SSDpositions
        if isfield(data.(file_list{x}), 'SSDpositions')
            for y = 1:number_pairs
                amplitude = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
                period    = data.(file_list{x}).SSD.combined_fits.period{1};
                phase     = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
                SSDlength_times  = data.(file_list{x}).SSDpositions.time_stamps;
                data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:) = Phases_SineFunc(SSDlength_times, period, amplitude, phase);
            end
        end
        for y = 1:number_fits
            if isfield(data.(file_list{x}), 'SSDpositions')
                %SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
                %data.(file_list{x}).SSD.surfmins_at_SSDlengths(y,:) = polyvalnm_solve2(surf_diff, 0, SSDlength_times);
                All_times_for_SSD_positions = data.(file_list{x}).SSDpositions.time_stamps - data.(file_list{x}).SSD.image_time(1,1);
                data.(file_list{x}).SSDpositions.surfmins_at_SSDlengths(y,:) = polyvalnm_solve2(surf_diff, 0, All_times_for_SSD_positions);
            end
        end
    end



    %get reference lengths and changes in length
    ref_positions = data.(file_list{x}).SSD.single_fits.reference_position{1};
    change_position = data.(file_list{x}).SSD.single_surfmins(:,1) - data.(file_list{x}).SSD.single_surfmins(:,end);
    
    for y = 1:number_pairs
        
        period_all(x,y) = data.(file_list{x}).SSD.combined_fits.period{1};
        period_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.period_err{1};
        
        phases_all(x,y) = data.(file_list{x}).SSD.combined_fits.phases{1}(y);
        phases_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.phases_err{1}(y);
        
        amplitudes_all(x,y) = data.(file_list{x}).SSD.combined_fits.amplitudes{1}(y);
        amplitudes_all_err(x,y) = data.(file_list{x}).SSD.combined_fits.amplitudes_err{1}(y);
        
        ref_length(x,y) = ref_positions(y+1) - ref_positions(y);

        change_length(x,y) = change_position(y+1) - change_position(y);
        if strcmpi(data.(file_list{x}).DispData.output_vars.length_type, 'pair') ==1
            change_length(x,y) = change_position(number_pairs+y) - change_position(y);
        elseif strcmpi(data.(file_list{x}).DispData.output_vars.length_type, 'anelastic') ==1
            change_length(x,y) = change_position(y+1) - change_position(y);
        else
            error('not implemented')
        end

        bulk_strain(x,y) = change_length(x,y)./ref_length(x,y);
        elapsed_time(x) = times2(end)-times2(1);
        bulk_strain_rate(x,y) = bulk_strain(x,y)./(times2(end)-times2(1));
        
        number_images(x) = data.(file_list{x}).used_images;
    end
         %keyboard   
end





%% plot stuff
disp('plotting stuff')

% 
% 
% %% plot position changes for each foil in turn.
% for x = 1 : size(file_list)    
%     if number_pairs <= 5
%                 
%         %sort data for plotting.
%         All_times      = data.(file_list{x}).image_time - data.(file_list{x}).image_time(1); 
%         All_positions  =  data.(file_list{x}).Disp.position_change + mean(data.(file_list{x}).reference_positions(1:2,:));
%         
%         %All_BackG   = data.(file_list{x}).Disp.BackGround; %in disp calculation bg is removed before sine is fit.
%         if isfield(data.(file_list{x}), 'Disp')
%             if isfield(data.(file_list{x}).Disp, 'sinusoid')
%                 Disp_times         = data.(file_list{x}).DispData.image_time;
%                 Disp_times_more    = data.(file_list{x}).DispData.times_more;
%                 Disp_sinusoid      = data.(file_list{x}).DispData.single_sinusoid;
%                 Disp_sinusoid_more = data.(file_list{x}).DispData.single_sinusoid_more;
% 
%                 Disp_ref_length    = data.(file_list{x}).DispData.output_values.reference_length{1};
% 
%                 Disp_backG         = data.(file_list{x}).DispData.BackGround'+Disp_ref_length';
%                 Disp_backG_more = [];
%                 for y = 1:size(Disp_backG,1)
%                     Disp_backG_more(y,:) = interp1(Disp_times, Disp_backG(y,:), data.(file_list{x}).Disp.times_more);
%                 end
% 
%                 Disp_lengths       = Disp_sinusoid +Disp_backG;
%                 Disp_lengths_more  = Disp_sinusoid_more+Disp_backG_more;
% 
%                 All_Fitted = ismember(All_times, Disp_times);
%            
%             end
%         else
%             All_Fitted = 1:size(data.(file_list{x}).image_time);
%         end
%         
%         if isfield(data.(file_list{x}), 'SSD') 
%             SSD_times              = data.(file_list{x}).SSD.image_time;
%             SSD_times_more         = data.(file_list{x}).SSD.times_more;%image_time;
%             SSD_sinusoid      = data.(file_list{x}).SSD.sinusoid;
%             SSD_sinusoid_more = data.(file_list{x}).SSD.sinusoid_more;%sinusoid;
% 
%             SSD_ref_length    = data.(file_list{x}).SSD.combined_fits.reference_length{1};
% 
%             if strcmpi(data.(file_list{x}).Disp.output_vars.length_type, 'pair')==1
%                 SSD_backG         = data.(file_list{x}).SSD.surfmins(1:end/2,:)-...
%                     data.(file_list{x}).SSD.surfmins(end/2+1:end,:)+ ...
%                     data.(file_list{x}).SSD.single_fits.reference_position{1}(end/2+1:end)-...
%                     data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end/2);
% 
%                 SSD_backG_more    = data.(file_list{x}).SSD.surfmins_more(1:end/2,:)-...
%                     data.(file_list{x}).SSD.surfmins_more(end/2+1:end,:)+...
%                     (data.(file_list{x}).SSD.single_fits.reference_position{1}(end/2+1:end))-...
%                     (data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end/2));
% 
%                 SSD_backG2 = CombineBoxes(data.(file_list{x}).Disp.output_vars.length_type, 'displacements', data.(file_list{x}).SSD.surfmins);
%                 keyboard
% 
%             elseif strcmpi(data.(file_list{x}).Disp.output_vars.length_type, 'anelastic')==1
%                 SSD_backG         = data.(file_list{x}).SSD.surfmins(1:end-1,:)-data.(file_list{x}).SSD.surfmins(2:end,:) + ...
%                                     data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end) -...
%                                     data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1);
%                 SSD_backG_more    = data.(file_list{x}).SSD.surfmins_more(1:end-1,:)-data.(file_list{x}).SSD.surfmins_more(2:end,:) + ...
%                                     data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end) - ...
%                                     data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1);
% 
%             else
%                 error('not implemented')
%             end
%             SSD_length        = SSD_sinusoid + SSD_backG;
%             SSD_length_more   = SSD_backG_more + SSD_sinusoid_more;
%             SSD_phases        = mod(data.(file_list{x}).SSD.phases + data.(file_list{x}).SSD.combined_fits.phases{1}', 2*pi);
%         end
%         
%         if isfield(data.(file_list{x}), 'SSDpositions')
%             SSDposition_backG         = data.(file_list{x}).SSDpositions.surfmins_at_SSDlengths(1:end-1,:)-data.(file_list{x}).SSDpositions.surfmins_at_SSDlengths(2:end,:)+...
%                                         (data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end))-(data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1));
%             SSDlengths_all         = data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(1:end-1,:)-data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(2:end,:);
%             SSDlengths_all         = data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths;
%             SSDlengths_all = SSDlengths_all+SSD_ref_length';
%         end
%         
% 
%         %plot all foil positions and fits together.
%         figure(1)
%         clf
% %         make_it_tight = true;
% %         subplot2 = @(m,n,p) subtightplot(m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
%         
%         for y = 1:number_pairs
%             ax(y) = subplot(number_pairs,1,y); cla, hold on,
%             
%             %plot displacment data
%             plot(All_times, All_length(:,y), 'b.') %SSD minimia
%             
%             %plot disp fit
%             if exist('Disp_times_more', 'var')
%                 plot(Disp_times_more, Disp_lengths_more(y,:), '-r') %displacement fit
%             end
%             if exist('SSD_times_more', 'var')
%                 plot(SSD_times_more, SSD_length_more(y,:), '-k') %SSD minima fit
%             end
%             plot(SSD_times_more, SSD_backG_more(y,:), '--k') %SSD length background
%             %plot(max(SSD_times_more)/2, SSD_ref_length(y), 'o') %SSD reference length
%             ylabel('Length (pixels)')
%             xlabel('Time (s)')
%             hold off
%             limits(:,y) = get(ax(y), 'Ylim');
%             
%             box on
%         end
%         
%         
%         %set y-axes to be a nice ratio of each other
%         if 1
%             extent = limits(2,:) - limits(1,:);
%             ratio = max(extent(:))./extent;
%             rr = floor(ratio);
%             rr_diffs = ratio-rr;
%             for y = 1:number_pairs
%                 if rr_diffs(y) < 0
%                     
%                     %mean 
%                     ax_mean = mean(limits(:,y));
%                     %range
%                     rng = max(extent(:));%./(rr(y)-1);
%                     
%                     ax_min = ax_mean-rng/2;
%                     ax_max = ax_mean+rng/2;
%                     
%                     set(ax(y), 'Ylim', [ax_min, ax_max]);
%                 elseif rr_diffs(y) > 0
%                     % expand the axes.
%                     
%                     %mean 
%                     ax_mean = mean(limits(:,y));
%                     %range
%                     rng = max(extent(:));%./rr(y);
%                     
%                     ax_min = ax_mean-rng/2;
%                     ax_max = ax_mean+rng/2;
%                     
%                     set(ax(y), 'Ylim', [ax_min, ax_max]);
%                 end
%             end
%         end
% 
%         %move figures to be close to each other. Remove space
%         pos=get(ax,'position');
%         bottom=pos{end}(2);
%         top=pos{1}(2)+pos{1}(4);
%         plotspace=top-bottom;
%         for y = length(pos):-1:1
%             pos{y}(2)=bottom+(length(pos)-y)*plotspace/length(pos);
%             pos{y}(2)=bottom+(length(pos)-y)*pos{y}(4);
%             %pos{y}(2) = bottom-pos{y}(4);
%             set(ax(y),'position',pos{y});
%         end
%         for y = 1:length(pos)
%             if rem(y,2) ~= 0
%                 set(ax(y),'YAxisLocation','right');
%             end
%             if y< length(pos)
%                set(ax(y),'xticklabel',[]);
%             end
%         end
%         
%         %label the axes
%         if 1
%             for y = 1:length(ax)
%                 if y==s
%                     title_str = sample_name;
%                 elseif y==e
%                     title_str = standard_name;
%                 else
%                     title_str = num2str(y);
%                 end
%                 tit(y)=title(ax(y),['  ',char(96+y),'. ',title_str]);
%                 tit(y).Position = [ ax(y).XLim(1), ax(y).YLim(2) ];
%                 tit(y).HorizontalAlignment = 'left';
%                 tit(y).VerticalAlignment = 'top';
%             end
%         end
%             
%         
%         
%         
%         %tightfig
%         saveimage(['Fits_',file_list{x}])
%         %keyboard
%         
%         %plot individual fits and residuals
%         % plot as tag group.
%         if 1
%             mainfig = figure(2);
%             tabgroup = uitabgroup(mainfig, 'Position', [.01 .01 .98 .98]);
%             
%             for y = 1:number_pairs
% 
%                 
%                 tab(y)=uitab(tabgroup,'Title', sprintf('Sample_%i', y));
%                 axes('parent',tab(y))
%                                 
%                 subplot(2,1,1), hold on,
% 
%                 %plot displacment data
%                 plot(All_times, All_length(:,y), 'b.')
%                 %plot disp fits
%                 if exist('Disp_times_more', 'var')
%                     plot(Disp_times_more, Disp_lengths_more(y,:), '-r') %displacement fit
%                 end
%                 if exist('SSD_times_more', 'var')
%                     plot(SSD_times_more, SSD_length_more(y,:), '-k') %SSD minima fit
%                 end
%                 
%                 if isfield(data.(file_list{x}), 'SSDpositions')
%                     SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
%                     windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
%                     plot(SSDlength_times', windowed_lengths, 'g.')
% 
%                     %SSD_all_ref = polyvalnm_solve2(squeeze(data.(file_list{x}).SSD.surf_diff(y,:,:)), 0, SSDlength_times);
%                     %SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);
%                     %plot(SSDlength_times', SSDposition_backG(y,:), 'g.-')
%                     %plot(SSDlength_times', SSD_all_ref, 'c.-')
%                 end
%                 hold off
%                 ylabel('Length change (pixels)')
%                 xlabel('Time (s)')
%                 title('Fits')
%                 x_lims = get(gca,'Xlim');
%                 
%                 subplot(2,1,2), hold on,
% 
%                 %plot residuals
%                 %plot(All_times(ismember(round(All_times,2),round(SSD_times(1,:),2))), All_length(ismember(round(All_times,2),round(SSD_times(1,:),2)),y)-SSD_length(y,:)', 'b.')
%                 plot(All_times, All_length(:,y)-SSD_length(y,:)', 'b.')
% 
%                 %plot disp fits
%                 if exist('Disp_times_more', 'var')
%                     plot(Disp_times_more, Disp_lengths_more(y,:)-SSD_length_more(y,:), '-r')
%                 end
%                 if exist('SSD_times_more', 'var')
%                     plot(SSD_times_more, SSD_length_more(y,:)-SSD_length_more(y,:), '-k')
%                 end
%                 
%                 if isfield(data.(file_list{x}), 'SSDpositions')
%                     SSDlength_times  = data.(file_list{x}).SSDpositions.time_stamps;
%                     windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
%                     %all fit at times of windowed data
%                     %background
%                     SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);
% 
%                     plot(SSDlength_times', windowed_lengths'-SSD_all_ref, 'g.')
% 
%                 end
%                 ylabel('Residuals (pixels)')
%                 xlabel('Time (s)')
%                 title('Residuals')
%                 set(gca, 'Xlim', x_lims);
%                 %set(gca, 'Ylim', [-0.05, 0.05])
% 
%     %             if y == 1
%     %                 title(file_list{x}, 'Interpreter', 'none')
%     %             end
%             end
%         end
% 
%         %% plot residuals for different fits (if applicable)
%         figure
%         for y = 1:number_pairs
%             ax(y) = subplot(1,number_pairs,y); cla, hold on,
% 
% 
%             plot(ax(y),[min(SSD_length(y,:)),max(SSD_length(y,:))], [min(SSD_length(y,:)),max(SSD_length(y,:))], 'k-')
% 
% 
%             plot(ax(y), SSD_length(y,:)', All_length(:,y), 'b.')
%             hold on
%             if isfield(data.(file_list{x}), 'SSDpositions')
%                 SSDlength_times  = data.(file_list{x}).SSDpositions.time_stamps;
%                 windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
%                 %all fit at times of windowed data
%                 %background
%                 SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);
% 
%                 plot(ax(y),SSD_all_ref', windowed_lengths', 'g.')
% 
%             end
% 
%             ylabel('point-by-point length')
%             xlabel('SSD lengths')
%             title('')
%         end
% 
% 
% 
%         figure
%         for y = 1:number_pairs
%             ax(y) = subplot(1,number_pairs,y); cla, hold on,
% 
%             resid = All_length(:,y)-SSD_length(y,:)';
%             disp('point-by point residuals')
%             disp(['mean ', num2str(mean(resid))])
%             disp(['std ', num2str(std(resid))])
%             plot(ax(y), sort(resid), (1:length(resid))/length(resid)*100, 'b.')
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
% 
% 
%             ylabel('Cumulative distribution (%)')
%             xlabel('Residual')
%             title('Residuals')
%         end
% 
% 
%     end
% end
% 



%% plot length changes
for x = 1 : size(file_list)    
    if 1;%number_pairs <= 4
        
        
        %sort data for plotting.
        All_times   = data.(file_list{x}).image_time - data.(file_list{x}).image_time(1); 
        All_length  = data.(file_list{x}).Disp.length_change + data.(file_list{x}).Disp.reference_length;
        
        %All_BackG   = data.(file_list{x}).Disp.BackGround; %in disp calculation bg is removed before sine is fit.
        if isfield(data.(file_list{x}), 'Disp')
            if isfield(data.(file_list{x}).Disp, 'sinusoid')
                Disp_times         = data.(file_list{x}).DispData.image_time;
                Disp_times_more    = data.(file_list{x}).Disp.times_more;
                Disp_sinusoid      = data.(file_list{x}).Disp.sinusoid;
                Disp_sinusoid_more = data.(file_list{x}).Disp.sinusoid_more;

                Disp_ref_length    = data.(file_list{x}).DispData.output_values.reference_length{1};

                Disp_backG         = data.(file_list{x}).DispData.BackGround'+Disp_ref_length';
                Disp_backG_more = [];
                for y = 1:size(Disp_backG,1)
                    Disp_backG_more(y,:) = interp1(Disp_times, Disp_backG(y,:), data.(file_list{x}).Disp.times_more);
                end

                Disp_lengths       = Disp_sinusoid +Disp_backG;
                Disp_lengths_more  = Disp_sinusoid_more+Disp_backG_more;

                All_Fitted = ismember(All_times, Disp_times);
           
            end
        else
            All_Fitted = 1:size(data.(file_list{x}).image_time);
        end
        
        if isfield(data.(file_list{x}), 'SSD') 
            SSD_times              = data.(file_list{x}).SSD.image_time;
            SSD_times_more         = data.(file_list{x}).SSD.times_more;%image_time;
            SSD_sinusoid      = data.(file_list{x}).SSD.combined_sinusoid;
            SSD_sinusoid_more = data.(file_list{x}).SSD.combined_sinusoid_more;%sinusoid;

            SSD_ref_length    = data.(file_list{x}).SSD.combined_fits.reference_length{1};

            if strcmpi(data.(file_list{x}).DispData.output_vars.length_type, 'pair')==1
                SSD_backG         = data.(file_list{x}).SSD.surfmins(1:end/2,:)-...
                    data.(file_list{x}).SSD.surfmins(end/2+1:end,:)+ ...
                    data.(file_list{x}).SSD.single_fits.reference_position{1}(end/2+1:end)-...
                    data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end/2);

                SSD_backG_more    = data.(file_list{x}).SSD.surfmins_more(1:end/2,:)-...
                    data.(file_list{x}).SSD.surfmins_more(end/2+1:end,:)+...
                    (data.(file_list{x}).SSD.single_fits.reference_position{1}(end/2+1:end))-...
                    (data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end/2));

                SSD_backG2 = CombineBoxes(data.(file_list{x}).Disp.output_vars.length_type, 'displacements', data.(file_list{x}).SSD.surfmins);
                keyboard

            elseif strcmpi(data.(file_list{x}).DispData.output_vars.length_type, 'anelastic')==1
                SSD_backG         = data.(file_list{x}).SSD.single_surfmins(1:end-1,:)-data.(file_list{x}).SSD.single_surfmins(2:end,:) + ...
                                    data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end) -...
                                    data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1);
                SSD_backG_more    = data.(file_list{x}).SSD.single_surfmins_more(1:end-1,:)-data.(file_list{x}).SSD.single_surfmins_more(2:end,:) + ...
                                    data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end) - ...
                                    data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1);

            else
                error('not implemented')
            end
            SSD_length        = SSD_sinusoid + SSD_backG;
            SSD_length_more   = SSD_backG_more + SSD_sinusoid_more;
            SSD_phases        = mod(data.(file_list{x}).SSD.combined_phases + data.(file_list{x}).SSD.combined_fits.phases{1}', 2*pi);
        end
        
        if isfield(data.(file_list{x}), 'SSDpositions')
            SSDposition_backG         = data.(file_list{x}).SSDpositions.surfmins_at_SSDlengths(1:end-1,:)-data.(file_list{x}).SSDpositions.surfmins_at_SSDlengths(2:end,:)+...
                                        (data.(file_list{x}).SSD.single_fits.reference_position{1}(2:end))-(data.(file_list{x}).SSD.single_fits.reference_position{1}(1:end-1));
            SSDlengths_all         = data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(1:end-1,:)-data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(2:end,:);
            SSDlengths_all         = data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths;
            SSDlengths_all = SSDlengths_all+SSD_ref_length';
        end
        
%         if isfield(data.(file_list{x}), 'SSDpositions')
%             keyboard
% 
%             data.(file_list{x}).SSDpositions.time_stamps
% 
%             SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
%             SSD_all_ref = polyvalnm_solve2(surf_diff, 0, SSDlength_times);
% 
%         end

        
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
            if exist('Disp_times_more', 'var')
                plot(Disp_times_more, Disp_lengths_more(y,:), '-r') %displacement fit
            end
            if exist('SSD_times_more', 'var')
                plot(SSD_times_more, SSD_length_more(y,:), '-k') %SSD minima fit
            end
            plot(SSD_times_more, SSD_backG_more(y,:), '--k') %SSD length background
            %plot(max(SSD_times_more)/2, SSD_ref_length(y), 'o') %SSD reference length
            ylabel('Length (pixels)')
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
            
        
        
        
        %tightfig
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

                    %SSD_all_ref = polyvalnm_solve2(squeeze(data.(file_list{x}).SSD.surf_diff(y,:,:)), 0, SSDlength_times);
                    %SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);
                    %plot(SSDlength_times', SSDposition_backG(y,:), 'g.-')
                    %plot(SSDlength_times', SSD_all_ref, 'c.-')
                end
                hold off
                ylabel('Length change (pixels)')
                xlabel('Time (s)')
                title('Fits')
                x_lims = get(gca,'Xlim');
                
                subplot(2,1,2), hold on,

                %plot residuals
                %plot(All_times(ismember(round(All_times,2),round(SSD_times(1,:),2))), All_length(ismember(round(All_times,2),round(SSD_times(1,:),2)),y)-SSD_length(y,:)', 'b.')
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
                %set(gca, 'Ylim', [-0.05, 0.05])

    %             if y == 1
    %                 title(file_list{x}, 'Interpreter', 'none')
    %             end
            end
        end

        %% plot residuals for different fits (if applicable)
        figure
        for y = 1:number_pairs
            ax(y) = subplot(1,number_pairs,y); cla, hold on,


            plot(ax(y),[min(SSD_length(y,:)),max(SSD_length(y,:))], [min(SSD_length(y,:)),max(SSD_length(y,:))], 'k-')


            plot(ax(y), SSD_length(y,:)', All_length(:,y), 'b.')
            hold on
            if isfield(data.(file_list{x}), 'SSDpositions')
                SSDlength_times  = data.(file_list{x}).SSDpositions.time_stamps;
                windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
                %all fit at times of windowed data
                %background
                SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);

                plot(ax(y),SSD_all_ref', windowed_lengths', 'g.')

            end

            ylabel('point-by-point length')
            xlabel('SSD lengths')
            title('')
        end



        figure
        for y = 1:number_pairs
            ax(y) = subplot(1,number_pairs,y); cla, hold on,

            resid = All_length(:,y)-SSD_length(y,:)';
            disp('point-by point residuals')
            disp(['mean ', num2str(mean(resid))])
            disp(['std ', num2str(std(resid))])
            plot(ax(y), sort(resid), (1:length(resid))/length(resid)*100, 'b.')
            hold on
            if isfield(data.(file_list{x}), 'SSDpositions')
                SSDlength_times              = data.(file_list{x}).SSDpositions.time_stamps;
                SSD_all_ref = SSDposition_backG(y,:)+data.(file_list{x}).SSDpositions.sinusoid_at_SSDlengths(y,:);
                %SSD_all_ref = SSD_all_ref+;
                windowed_lengths = data.(file_list{x}).SSDlengths.foilpositions(:,y);
                resid_windowed = windowed_lengths'-SSD_all_ref;%SSD_all_ref-windowed_lengths'+mean(windowed_lengths);


                plot(ax(y), sort(resid_windowed), (1:length(resid_windowed))/length(resid_windowed)*100, 'g.')


                disp('windowed residuals')
                disp(['mean ', num2str(mean(resid_windowed))])
                disp(['std ', num2str(std(resid_windowed))])


            end


            ylabel('Cumulative distribution (%)')
            xlabel('Residual')
            title('Residuals')
        end





        %% plot stress--strain
        if number_fits ==3
            
            
            figure(3); clf
            hold on,
            colormap('hsv')
            plot((All_length(All_Fitted,s)-SSD_backG(s,:)')/SSD_ref_length(s), (All_length(All_Fitted,e)-SSD_backG(e,:)')/SSD_ref_length(e), ':k')
            scatter((All_length(All_Fitted,s)-SSD_backG(s,:)')/SSD_ref_length(s), (All_length(All_Fitted,e)-SSD_backG(e,:)')/SSD_ref_length(e), 9, SSD_phases(s,:), 'filled');
            if exist('Disp_times_more', 'var')
                plot(Disp_sinusoid_more(s,:)/Disp_ref_length(s), Disp_sinusoid_more(e,:)/Disp_ref_length(e), '-r', 'LineWidth', 1)
            end
            if exist('SSD_times_more', 'var')
                plot(SSD_sinusoid_more(s,:)/SSD_ref_length(s), SSD_sinusoid_more(e,:)/SSD_ref_length(e), '-k', 'LineWidth', 1)
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
            %tightfig
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
    
    
