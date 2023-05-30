%% Foil Track, Example 1.
%
% Running these steps should reproduce the files in the 'Outputs' directory. 
% 
% N.B. the values phase and amplitude values calculated by this script are not identical to those
% reported in the paper describing the method. 
% This is because the original images were in a NetCDF file, but to reduce the size of the example data
% the images were extracted and reduced in size. The time stamps for the extracted images have a resolution of 
% 1s (the time stamp of when the image was created/edited). This is less than the <0.01s resoultion for the time stamps
% in the original file. 
% In the example dataset the images have idea time stamps 2.5s apart, which cannot be maintained. So the times spacing 
% of the images is different. 
% 
% The exmaple outputs are calculated from the extracted example data and so should be reporducible. 
%
% Tested with Matlab 2021b on MacOS 16.6
%
% Simon Hunt, 2023



%% STEP 1: make the list file from the tif images. 
% To make a list file for the directory of image files:
MakeSingleLstFile('./data/Zn_08_27tons_117C_100s_027')

% For multiple directoryies call: 
% MakeManyTimes('./data/')
% both of these commands will make the file 'Zn_08_27tons_117C_100s_027.lst'
% which is a list of all the image files in the data directory.



%% STEP 2: Set up the experiment.
% To set the options for the analysis run:
AnalysisOptions

% The details of what each line sets are in the documentation and/or AnalysisOtions header

% To reproduce the output files in the example, the settings are: 
% Experiment name = 'Zn_08' -- this corresponds to the initial string of the file name and occurs in all the files in the experiment.
% experiment location = 'NSLS: X17B2' 
% reference images = "ref"
% Length type = "anelastic"
% all other settings are the defaults

% This makes the file 'Zn08.mat'



%% STEP 3: Select Boxes for analysis
% To select the regions of interest for the image processing:
ImageAnalysis('Zn_08_27tons_117C_100s_027.lst','boxes')
% this works through a number of steps. 
% a. The image rotation is 270 degrees
% b. To select individual boxes click top left and bottom right of each foil.
% c. work way through dialogues. 

% This script eventualy makes the file 'Zn_08_27tons_117C_100s_027_boxes.mat' and the image 'Zn_08_27tons_117C_100s_027_boxes.tif'

% The ultiamate output of the script depends on the postion of the boxes that are selected. 
% To make idential output files to those in the 'outputs' directory copy the '*_boxes.mat' file to the directory this file is in. 



%% STEP 4: Perform the image analysis
% To process all the images:
ImageAnalysis('Zn_08_27tons_117C_100s_027.lst','disp')

% this will produce a number of files:
% - 'Zn_08_27tons_117C_100s_027_position_change.txt' -- displacments of foils realtive to refernece, pair by pair.
% - 'Zn_08_27tons_117C_100s_027_SSD.mat' -- Matlab files with all the SSD displacements in. 



%% STEP 5a: Calculate phases of the deformation
% This can be done from the displacement data: 
Phases('Zn_08_27tons_117C_100s_027_position_change.txt', 'period', 100);

% or from the SSD met file: 
PhasesSSD('Zn_08_27tons_117C_100s_027_SSD.mat', 'period', 'title');

% for large data sets with multiple data at each period it is possible to check that all the periods are consistent:
PhaseCheck('type', 'disp')
PhaseCheck('type', 'ssd')

% this will produce a number of files:
% - 'Zn_08_27tons_117C_100s_027_sine_fits.txt' -- phase and ampltuide of the sample length change as defined in the ZN08.mnat file.
% - 'Zn_08_27tons_117C_100s_027_sine_fits_single.txt' -- phase and ampltuide for the fits of each foil separately.
% - 'Zn_08_27tons_117C_100s_027_sine_fits_FITS.mat' -- The values stored as matlab data file.
% - 'Zn_08_27tons_117C_100s_027_SSD_sine_fits.txt' -- The values stored as matlab data file.
% - 'Zn_08_27tons_117C_100s_027_SSD_FITS.mat' -- The phase and amplitude for each region of interest stored as matlab data file.



%% STEP 5b: Calculate positions and rates
% These scripds calculate the postion and position change-rates for the two
% data methods. 

% set the length of the window to work over (number of images)
window = 7;
%set the degree of the polynomial
degree = 2;

%run the script:
PositionTime('Zn_08_27tons_117C_100s_027_position_change.txt',  'window', window, 'overwrite', 'degree', degree)
% this will produce a number of files:
% - Zn_08_27tons_117C_100s_027_window_PositionChangeRate.txt
% - Zn_08_27tons_117C_100s_027_window_LengthChangeRate.txt

%run the scripts
PositionTimeSSD('Zn_08_27tons_117C_100s_027_SSD',  'window', window, 'surface coefficients', 6,degree, 'overwrite')
% this will produce a number of files:
% - Zn_08_27tons_117C_100s_027_SSDwindow_positions.txt
% - Zn_08_27tons_117C_100s_027_SSDwindow_lengths.txt
% - Zn_08_27tons_117C_100s_027_SSDwindow_PositionChangeRate.txt
% - Zn_08_27tons_117C_100s_027_SSDwindow_LengthChangeRate.txt




% =================
% Plot the outputs of the fitting. 
% =================

%% STEP 6a: phases 
% This script compares the phases and amplitudes of the two data fitting
% methods.
SSD_displacement_compare('Zn_08_27tons_117C_100s_027', 'samples', {'corundum, top','zinc','corundum, bottom'} )

% this will produce the remaining image files in the outputs directory. 
% these figures are equivalent to the figures in the manuscript describing this method. 
% - Fits_Zn_08_27tons_117C_100s_027.pdf
% - PhasePhaseZn_08_27tons_117C_100s_027.pdf
% - AmplitudeAmplitudeZn_08_27tons_117C_100s_027.pdf
% - PeriodPeriodZn_08_27tons_117C_100s_027.pdf
% - PhaseLagPhaseLagZn_08_27tons_117C_100s_027.pdf



%% STEP 6b: Calculate positions and rates
% This script compares the positions and position chagne rates of the two data fitting
% methods.
SSD_strainrate_compare('Zn_08_27tons_117C_100s_027', 'samples', {'corundum, top','zinc','corundum, bottom'} )

% this will produce the remaining image files in the outputs directory. 
% - StrainRateFits_Zn_08_27tons_117C_100s_027.pdf
% - RateResiduals_Zn_08_27tons_117C_100s_027.pdf
% - RateResiduals2_Zn_08_27tons_117C_100s_027.pdf
