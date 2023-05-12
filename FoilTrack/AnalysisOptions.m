function out = AnalysisOptions(varargin)
%AnalysisOptions
% Sets-up the options to analyse X-radiography images, for rheology,
% thermal diffusivity or anelasticity experiments.
%
% Green boxes are likely to need to be changed, orange possibly and red unlikely to need 
% changing.
%
% The options presented are:
%   Headers:-
%       Experiment name:    
%           This is the name of the experiment and the name of the saved options file. 
%           It should be the name of the experiment and the leading string in the image 
%           file names.
%       Experiment location: 
%           Pull down menu to select where the experiment was performed. This is needed
%           so the correct image preprocessing (rotation, inversions) can be performed. 
%           The 'sides' option (where present) is a 90 degree rotation of the images so 
%           that displcements of the sides of the foils can be determined. 
%
%   ImageAnalysis options. These are the options for the ImageAnalysis script that 
%           processes the images:-
%       Reference Image:    
%           'Sequence' uses each image and compares it to the following image, moving 
%           the reference box as it goes. 'Reference' compares the middle image in the sequence 
%           to all the other images. 
%       Search distance:
%           This is the distance (in pixels) that the code calculates the SSD either side 
%           of the reference position. The default for the X17B2 / 6-BM-B beamlines is 5.
%       Minimum SSD:
%           The method by which the minimum in the SSD is interpolated. The options are 
%           'spline', 'ploynomial', and 'gaussian'.
%       Output data:
%           What files the ImageAnalysis script saves. The 'displacements' is the interpolated 
%           displaceemnts, whilst SSD saves all the calcualted SSD in a single file. 
%       Bright spots:
%           If 'yes', this removes the very bright spots from the images. This is used 
%           when the images are all black but for the spots so the data cannot be seen. 
%           The option was use for the Princeton camera at the X17B2/X17B2ss beamlines 
%           which produced 8 bit data within a 12 or 16 bit file.
%       Mask background:
%           Tries to recognise the regions outside the foil but within the slected box 
%           and replaces them with NaN. This removes the effect of these pixels from the 
%           analysis. 
%
%   Phases options:-
%       Length type:
%           Determines how the script combines multiple foils. 'Single' treats each box individually,
%           'Pair' aligns boxes with the same horizontal position, with more than 2 foils it 
%           combines the two outside foils works inwards in pairs. 
%       Image Scaling:
%           Size of each pixel in microns. FIX ME. IS THIS CORRECT?
%       Discard:
%           If there length type is pair and there are an odd number of foils which one to discard. 
%       Symmetry type:
%           How to combine the foils for the pair length type. Is the
%           bottom foil switched round?
%       Background type:
%           How to detrend the data before fitting the sinusoid.
%       Errors:
%           Save these or not?

% Simon Hunt 2011 - 2016
% $ Version 1.1.3 $ 

% THIS NEEDS MANY THINGS DOING...
% - rheology/t.d. defaults adding as switch?
% - Kappa_solve settings?

% array of contents for 'figure' - allows for loops in creating options choices
beamlines = {...'PRESENT BUT NOT IMPLEMENTED';
    '';
    'APS: 6-BM-B';
    'APS: 13-ID-D';
    'APS: 13-ID-D sides';
    'APS: 16-ID-B';
    'APS: 16-ID-B xrd';
    'DLS: I12';
    'ESRF: ID06LVP';
    'DESY: MAX200x';
    'NSLS: X17B2';
    'NSLS: X17B2ss';
    'NSLS2: MaxPD';
    'generic';};
beamlines(:,2) = CleanBeamlineName(beamlines(:,1));
expt_loc = 1;

% array of variables to be set for X-radiography analysis via. ImageAnalysis and Phases.
% columns are: associatio, variable name, title for selection option, text descripton , options, default setting, background colour
UI_content = {};
UI_content(size(UI_content,1)+1,:) = {'disp',  'analysis_type', 'Reference image',                 'What sort of analysis is wanted?',                                     {'Seq', 'Ref'},                                 'Seq',          'green', 1};
% UI_content(size(UI_content,1)+1,:) = {'disp',  'box_fit_type','Select Boxes: Box Positions',      'What feature/method to is key in selecting the boxes',                 {'Min', 'Slope', 'Gauss', 'Fit','None'},       'min',          'green'};
UI_content(size(UI_content,1)+1,:) = {'disp',  'search',       'Displacement: How far to look?',   'How far from reference position should SSD be calculated (+/-)?',      {'3', '4', '5', '7', '10', '20'},                 '5',            'orange', 2};
UI_content(size(UI_content,1)+1,:) = {'disp',  'min_type',     'Displacement: Minimum of SSD',     'What method for finding minimum of SSD hence displacement',            {'Spline', 'Poly', 'Gauss'},                    'Poly',         'orange', 1};
UI_content(size(UI_content,1)+1,:) = {'disp',  'out_data_type','Displacement: Output data types',  'Output data type from displacement calculations',                      {'Displcements', 'SSD', 'Both'},                'Both',         'orange', 1};
UI_content(size(UI_content,1)+1,:) = {'disp',  'spot_removal', 'Displacement: Bright spots',       'Mask out bright spots from seleced areas?',                            {'Yes', 'No'},                                  'Yes',          'red', 1};
UI_content(size(UI_content,1)+1,:) = {'disp',  'NaN_bg',       'Displacement: Mask background?',   'Replace area not in foil with ''NaN'' (use for poor quality data)',    {'Yes', 'No'},                                  'No',           'red', 1};
UI_content(size(UI_content,1)+1,:) = {'disp',  'offset',       'Displacement: Whole image movement?',   'Find displacement of areas of interest in whole image first?',    {'Yes', 'No'},                                  'No',           'red', 1};

UI_content(size(UI_content,1)+1,:) = {'phase', 'length_type', 'Length type',        'Calculate the phases from pairs or single foils',                      {'Deformation', 'Pair', 'Single', 'Anelastic'},                'Pair',         'green', 1};
UI_content(size(UI_content,1)+1,:) = {'phase', 'scaling',     'Image scaling',      'Number pixels per micron',                                             {'1', '2', '3', '4', '5', '6'},                 '2',            'green', 3};
UI_content(size(UI_content,1)+1,:) = {'phase', 'excessive',   'Excessive',          'Ignore points more then n sigma from the fit',                         {'2', '7'},                                     '4',            'orange', 2};
UI_content(size(UI_content,1)+1,:) = {'phase', 'get_rid',     'Discard',            'In the case of odd numbers of foils select which one to remove',       {'End', 'Middle', '1', '2', '3', '4', '5'},     'Middle',       'red', 1};
UI_content(size(UI_content,1)+1,:) = {'phase', 'symm_type',   'Symmetry type',      'Length changes using symmetic or antisymmetric boxes',                 {'Symmetric', 'Antisymmetric'},                 'Symmetric',    'red', 1};
UI_content(size(UI_content,1)+1,:) = {'phase', 'bg_type',     'Background type',    'Choose method for background removal',                                 {'Moving Average', 'Spline', 'Polynomial'},     'Spline',       'red', 1};
UI_content(size(UI_content,1)+1,:) = {'phase', 'err_calc',    'Errors?',            'Calcualte the errors to the fits?',                                    {'Yes', 'No'},                                  'Yes',          'red', 1};


%generate a list of variable names
var_name_array = UI_content(:,2);


%define default settings array
for x = 1:length(var_name_array)
    var = UI_content{x,6};
    if isnumeric(var) == 0
        param.(var_name_array{x}) = var;
    else
        param.(var_name_array{x}) = num2str(var);
    end
end


%% check for and read settings .mat file.
if nargin >=1  %if the settings file is passed into AnalysisOptions.
    [~,varfilename,~] = fileparts(varargin{1});
else
    nc_files = dir('*.nc');
    mat_files = dir('*.mat');
    files = dir('*.*');
    if isempty(nc_files) == 0
        nc_files(1).name;
        out = FileTitleInformation(nc_files(1).name);
        varfilename = out{1};
    elseif isempty(mat_files) == 0 %FIX ME. If not the guess the name of the directory.
        for x = 1 : length(mat_files)
            v = whos('-file', mat_files(x).name);
            if any(strcmpi({v.name}, 'min_type')) %pick random variable that should be present
                varfilename = mat_files(x).name(1:end-4);
                %         try
                %             out = FileTitleInformation(mat_files(x).name);
                %             varfilename = out{1};
            end
        end
    else
        varfilename = '';
    end
end

% if isempty(varfilename) == 1
%     varfilename = '';
% end

% if the varible file exists then read it.
if exist([varfilename,'.mat'],'file') == 2
    
    disp([varfilename,'.mat already exists. Loading contents.']);
    
    loaded_param = load([varfilename,'.mat'], var_name_array{:});
    %format the values 
    for x = 1: length(var_name_array)
        if isfield(loaded_param, var_name_array{x})
            if isnumeric(loaded_param.(var_name_array{x})) == 0
                param.(var_name_array{x}) = loaded_param.(var_name_array{x});
            else
                param.(var_name_array{x}) = num2str(loaded_param.(var_name_array{x}));
            end
        end
    end 
    
    %load the experiment location.
    locat = load([varfilename,'.mat'], 'expt_location');
    if isfield(locat, 'expt_location')== 1
%         locat.expt_location = str2loc(locat.expt_location, beamlines); %add spaces colons and hyphen into the saved strings.
%         expt_location = find(~cellfun(@isempty,strfind(beamlines, locat.expt_location)));
        expt_loc = FindBeamline(beamlines, locat.expt_location);
    else
        expt_loc = [];
    end
    if isempty(expt_loc)
        expt_loc = 1;
    end
end


%% create figure
% set basic dimensions for all parts of the user interface. 

% number of rows and columns needed
names = unique(UI_content(:,1));
number_columns = length(names);
for x = 1 : number_columns
    number_rows(x) = sum(strcmp(UI_content(:,1),names(x)));
end
number_rows = max(number_rows);


% part dimensions
height = 15;
left_border = 10;
bottom_border = 5;
text_column = 350;
column = 50;
panel_border = .0125;
entries = (1 - 2*panel_border)/number_rows;

save_button_height = 30;
save_button_width = 120;
save_height = save_button_height + bottom_border;
number_save_buttons = 2;


% size of descriptive text for each option.-- these sizes go within button groups
button_description_size = [left_border bottom_border+height text_column height];

%size and location of the button group panels which select options -- these fit within pannels 
buttongroup_size = [   left_border    bottom_border   text_column+2*left_border  3*height+bottom_border ];
buttongroup_size = repmat(buttongroup_size,number_rows,1);
steps = buttongroup_size(1,4) + bottom_border/2;
buttongroup_size(:,2) = buttongroup_size(:,2) + steps * (number_rows-1 : -1 : 0)';

%panel sizes -- fits within figure and defines width of the figure.
pannel_size = [ left_border   bottom_border + save_height  buttongroup_size(1,3)+2*left_border  steps*number_rows+bottom_border+height];
pannel_size = repmat(pannel_size, number_columns,1);
pannel_size(:,1) = pannel_size(:,1) + (0: number_columns-1)' * (pannel_size(1,3) + pannel_size(1,1));

%top bar size 
%top_size(3) is the number which defines the width of the figure
top_size = [0  pannel_size(1,2)+pannel_size(1,4)  pannel_size(number_columns,1)+pannel_size(number_columns,3)+left_border 4*height ];

% OK, cancel button sizes and locations...
bottom_button_sizes = [ 0 bottom_border save_button_width save_button_height];
bottom_button_sizes = repmat(bottom_button_sizes, number_save_buttons,1);
bottom_button_sizes(:,1) = -(number_save_buttons-1 : -1 :0)' * (save_button_width + left_border);
bottom_button_sizes(:,1) = bottom_button_sizes(:,1) + repmat(top_size(3)-left_border-save_button_width,number_save_buttons,1);

% figure_size 
figure_position = [250 238 top_size(3) top_size(2)+top_size(4)];


% colours
myred = [225/255 183/255 183/255];
mygreen = [179/255 209/255 189/255];
myorange = [239/255 190/255   117/255];

            
poss = zeros(1,number_columns);     


            
% Build UI -- figure
if findobj('Type','figure','Tag','Settings Window') == 1
    fig = findobj('Type','figure','Tag','Settings Window');
else
    fig = figure('Name','X-radiography Analysis Options','NumberTitle','off', 'MenuBar', 'none', 'Tag', 'Settings Window');
end
dims = get(fig,'Position');
set(fig,'Position', figure_position, 'Resize', 'off');
dims = get(fig,'Position');

%top panel for file name and location 
pan_top  = uipanel('Parent',fig,'Units','Pixels', 'Position',top_size,'bordertype','none');
info1 =  uicontrol('Parent',pan_top, 'Style','text','String','Experiment name:','Position',[left_border 1.5*bottom_border+2*height text_column height],'HorizontalAlignment','left');
inform = uicontrol('Parent',pan_top, 'Style','edit','String',varfilename,'pos',[left_border+2*column 1.5*bottom_border+2*height 4*column height+3],'background','white', 'Callback',@edittext_name);
info2 =  uicontrol('Parent',pan_top, 'Style','text','String','Experiment location:','Position',[left_border 2*bottom_border text_column height],'HorizontalAlignment','left');
inform2 = uicontrol('Parent',pan_top, 'Style','popupmenu','String', beamlines(:,1),'pos',[left_border+2*column 2*bottom_border 4*column height+3],'background','white', 'Callback',@editlocation);
set(inform2, 'Value', expt_loc);

%pannels for options 
% FIX ME: this is a bit of a cheat by naming the pannels and defining variables here. Needs to be done at the top really.
pan_disp = []; pan_phase = []; pan_phase2 = [];
pan_names = {'ImageAnalysis options'; 'Phases options'; 'extra'};
for x = 1 : number_columns
    eval(['pan_',names{x,1},' = uipanel(''Parent'',fig,''Title'',pan_names{x,:},''Units'',''Pixels'', ''Position'',pannel_size(x,:));']);
end

bg = 1;
% create buttongroups
for x = 1: size(UI_content,1)
    %count to keep spacings for options
    opt = strcmp(UI_content(x,1),names(:))'; % which row current options go in in column
    poss = opt+poss; %keeps track of how many options in each column
    row = poss(find(opt==1)); %which row current option goes in

    % draw the button group itself.
    if UI_content{x,8} == 1 %button group options
        
        g(x) = uibuttongroup('Parent',eval(['pan_',UI_content{x,1}]), 'Units','Pixels', 'Position',buttongroup_size(row,:),...
            'Title', UI_content(x,3), 'Tag', UI_content{x,2}, 'Background', eval(strcat(['my',UI_content{x,7}])));
        info1 = uicontrol('Style','text','String',UI_content(x,4),'Units','Pixels','Position',button_description_size,'Parent',g(x),'HorizontalAlignment','left', 'Background', eval(strcat(['my',UI_content{x,7}])));
        % Create radio buttons in the button group.

        for y = 1 : length(UI_content{x,5}) 
            vals = UI_content{x,5};
            u(x,y) = uicontrol('Style','Radio','String',vals(y),'pos',[left_border+column*(y-1) bottom_border column height],'Parent',g(x),'HandleVisibility','off', 'Background', eval(strcat(['my',UI_content{x,7}])));
        end
        
        % Initialize some button group properties.    
        def = strcmpi(vals, param.(UI_content{x,2}));
        def = (def == 1);
        set(g(x),'SelectionChangeFcn',@selcbk);
        set(g(x),'SelectedObject',u(x,def));  % Set defaults or previous values
        
    elseif UI_content{x,8} == 2 %slider 
        vals = UI_content{x,5};
        
        %panel
        pan_mid(x)  = uipanel('Parent',eval(['pan_',UI_content{x,1}]),'Units','Pixels', 'Position',buttongroup_size(row,:),...   );%,'bordertype','none');
                                'Title',UI_content(x,3),'Background', eval(strcat(['my',UI_content{x,7}])));
                            
        %text above slider
        info1 = uicontrol('Style','text','String',UI_content(x,4),'Units','Pixels','Position',button_description_size,'Parent',pan_mid(x),...
            'HorizontalAlignment','left', 'Background', eval(strcat(['my',UI_content{x,7}])));
        
        %slider
        steps = 1/(str2num(vals{end}) - str2num(vals{1}));
        g(x) = uicontrol('Style','Slider','pos',[left_border+column/2 bottom_border column*4 height],...
            'Min', str2num(vals{1}), 'Max', str2num(vals{end}), 'Value', str2num(param.(UI_content{x,2})), 'Tag', UI_content{x,2},...
            'Parent',pan_mid(x),'HandleVisibility','off', 'Background', 'w', 'Callback', @sliderbk, 'SliderStep',[steps steps]); 
        
        %text to right of slider
        val_unfo.(UI_content{x,2}) = uicontrol('Style','text','String',param.(UI_content{x,2}),'Units','Pixels','Position',[column*5 bottom_border column height] ,'Parent',pan_mid(x),...
            'HorizontalAlignment','left', 'Background', eval(strcat(['my',UI_content{x,7}])));      

    elseif UI_content{x,8} == 3 %string entry
        
        %panel
        pan_mid(x)  = uipanel('Parent',eval(['pan_',UI_content{x,1}]),'Units','Pixels', 'Position',buttongroup_size(row,:),... %,'bordertype','none');
                                'Title',UI_content(x,3),'Background', eval(strcat(['my',UI_content{x,7}])));
        
        %text above string entry 
        info1 = uicontrol('Style','text','String',UI_content(x,4),'Units','Pixels','Position',button_description_size,'Parent',pan_mid(x),...
            'HorizontalAlignment','left', 'Background', eval(strcat(['my',UI_content{x,7}])));
        
        %string entry
%         info1 =  uicontrol('Parent',pan_top, 'Style','text','String','Experiment name:','Position',[left_border 1.5*bottom_border+2*height text_column height],'HorizontalAlignment','left');
        g(x) = uicontrol('Parent',pan_mid(x), 'Style','edit','String',num2str(param.(UI_content{x,2})),'Tag',UI_content{x,2},...
                            'pos',[left_border+column/2 bottom_border column*4 height],'background','white',...[left_border 1.5*bottom_border+2*height column height+3],'background','white',...
                            'HorizontalAlignment','left', 'Callback',@textbk);
    else
        error('This call back type option is not recognised.');
                        
    end
    
end


% Create OK & Cancel push buttons.
pbh1 = uicontrol('Style','pushbutton','String','OK','Units','Pixels','Position',bottom_button_sizes(1,:),'Callback',@OKButtonCallback,'Parent',fig);
pbh2 = uicontrol(fig,'Style','pushbutton','String','Cancel','Units','pixels','Position',bottom_button_sizes(2,:),'Callback',@CancelButtonCallback,'Parent',fig);
%Create manytimes files buttons.
pbh3 = uicontrol('Style','pushbutton','String','Create batch files','Units','normalized','Position',[.0125 .0125 .2 .06],'Callback',@ManytimesButtonCallback,'Parent',fig);

%stops fending script until #ok or #cancel are pressed.
uiwait
out = [];
pause(1)

%% Callback functions

    %change variables for button groups
    function selcbk(source,eventdata) %#ok<INUSD>
        varnam = genvarname(get(source,'Tag'));
        param.(varnam) = char(get(eventdata.NewValue,'String'))
%         eval([varnam,' = get(eventdata.NewValue,''String'')'])
    end

    function sliderbk(hObj,~,~) %#ok<INUSL> %change variable for sliders
    % Called when user moves the slider control
    varnam = genvarname(get(hObj,'Tag'));
    varval = round(get(hObj,'Value'));
    param.(varnam) = num2str(varval)
    
    %set the string and slider to the rounded value.
    set(val_unfo.(varnam), 'String', varval);
    set(hObj, 'Value', varval);
    
    end
    
    
    function textbk(source, eventdata)%  change variable for strings
        varnam = get(source,'Tag');
        varval = get(source,'String');
        param.(varnam) = num2str(varval)
    end
  

    % location 
    function editlocation(source, ~)
%         val = get(source,'Value');
        expt_loc = get(source,'Value');
        expt_location = beamlines(expt_loc,2)
    end

    % Cancel button
    function CancelButtonCallback(~, ~)
        uiresume
        delete(fig);
    end



    % OK button
    function OKButtonCallback(~, ~)
        if isempty(varfilename) == 1
            warndlg('The experiment name is not set');
        else
            uiresume;
            delete(fig);
            
            % sets options to correct format. If all possibilities for options are numeric leave alone otherwise make lowercase character string
            fields = fieldnames(param);
            for count = 1:length(fieldnames(param))
                if isempty( str2num(param.(fields{count})) ) == 1
                    param.(fields{count}) = lower(param.(fields{count}));
                else
                    param.(fields{count}) = str2num(param.(fields{count}));
                end
            end
            
            expt_location = beamlines{expt_loc,2};
                        
            %saves the output file.
            if exist([varfilename,'.mat'],'file') == 0
                save(varfilename, '-struct', 'param');
                save(varfilename, 'expt_location', '-append');
            else
                save(varfilename, '-struct', 'param', '-append');
                save(varfilename, 'expt_location', '-append');
            end          
           
        end
    end


    % Mantyimes button
    function ManytimesButtonCallback(~, ~)
        makemanytimes
        
    end

    % file name 
    function edittext_name(source, ~)
        varfilename = get(source,'String')
    end

end


function expt_location = CleanBeamlineName(expt_location)

%remove the spaces, hyphens and colons from the experiment location
%strings.
expt_location = strrep(expt_location, ' ', '');
expt_location = strrep(expt_location, ':', '_');
expt_location = strrep(expt_location, '-', '');
end


function expt_location = FindBeamline(beamlines, loc)

expt_location = find(strcmpi(beamlines(:,2),loc)==1);

end

% 
% function names = loc2str(names)
% 
% %remove the spaces, hyphens and colons from the experiment location
% %strings.
% for x = 1:numel(names)
%     names{x} = strrep(names{x}, ' ', '');
%     names{x} = strrep(names{x}, ':', '_');
%     names{x} = strrep(names{x}, '-', '');
% end
% end
% 
% function location = str2loc(location, possibilities)
% 
% if ischar(location)
%     
%     %assumes that the saved location exists on the list of possibilities
%     possible = loc2str(possibilities);
%     location = loc2str(location);
%     try
%         location = possibilities{~cellfun(@isempty,strfind(possible,location))};
%     catch
%         location = possibilities{1};
%     end
% else
%     location = possibilities{location};
% end
%     
% end

% Versions
%  - 1.1.3 - 20th June 2016
%       - Fixed the experiment location so that it does not save a number as the location. 
%  - 1.1.2 - 3nd June 2016
%       - Added input options so that can pass the options file name into this script.
%  - 1.1.1 - 2nd June 2016
%       - Edits to mat file read to account for the fact that FileTitleInformation no longer throws out errors.
%  - 1.1 - 18th May 2016
%       - Locaitons implemented in the code. Change the string structioned
%       saved so that it can be readily recognised by ImageAnalysis.
%  - 1.0.1 - 7th May 2016
%       - edit to stop crashing if filename does not exist. (line 134).
%  - 1.0 - April 2016
%       - rewritten so variables are not unique variables but part of
%       structure. Therefore do not have to be predefined. 
%  - 0.5 - Jan 2014
%       - rewritten so nunber of options determines the size of the figure.
%  - 0.4 - Jan 2014
%       - generalised the code so can add new variables by adding lines to initialisation arrays.
%  - 0.3 - January 2012
%       - options presented are now read from 'expt'.mat if it exists.
%  - 0.2a
%		- added spot_removal to options
%  - 0.2 
%		- added implementation for options for ImageAnalysis.m
%  - 0.1 
%		- implemented options for Phases.m