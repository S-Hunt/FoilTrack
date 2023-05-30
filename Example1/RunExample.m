


%% make the list file from the tif images. 

MakeManyTimes


MakeSingleLstFile('../Zn_08_27tons_117C_100s_027')

% this command makes the file Zn_08_27tons_117C_100s_027.lst
% this file lists all the fit files in the data directory. 


%% Set up the experiment.

AnalysisOptions

% for more details on how to set up the expeirment see ... docuents. 

% here we need the options:
%Experiment name = 'Zn_08' -- this corresponds to the initial string of the
%file name and occurs in all the files in the experiment.

% experiment location = 'NSLS: X17B2' 

% reference images = "ref"

%Phase options // length type = "anelastic"

%everthing else is the default.

%clicking OK makes file Zn08.mat

%% Select Boxes for analysis

ImageAnalysis('Zn_08_27tons_117C_100s_027.lst','boxes')


%rotate 270 degrees

%select boxes. click top left and bottom right of each foil.
% work way through pictures

% makes boxes file.


%% do boxes

ImageAnalysis('Zn_08_27tons_117C_100s_027.lst','disp')



%% phases

Phases('Zn_08_27tons_117C_100s_027_position_change.txt', 'period', 100);

PhaseCheck('type', 'disp')


PhasesSSD('Zn_08_27tons_117C_100s_027_SSD.mat', 'period', 'title');

PhaseCheck('type', 'ssd')



%% plot outputs.

SSD_displacement_compare('Zn_08_27tons_117C_100s_027', 'samples', {'corundum, top','zinc','corundum, bottom'} )