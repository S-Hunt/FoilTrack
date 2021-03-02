
function choice = SelectBoxes_MultiButton(title, message, options)

% title = 'Hello';
% message = 'Does this bloody work?';
% options = {'Yes', 'No', 'Yeah right!','My Arse it does'};

%define some lengths
box_width = 200;
button_size = [100 30];
button_step = button_size(2)*1.2;

% Define width of text relative to box width.
colw = box_width-30;

%make diagloue box
FigPos    = get(0,'DefaultFigurePosition');
FigPos(3) = box_width;
FigPos(4) = 100;
FigPos    = getnicedialoglocation(FigPos, get(0,'DefaultFigureUnits'));
% set(d,'Position',FigPos)

% d = figure('Name',title,'Position',FigPos);
% d = dialog('Name',title,'Position',FigPos);
d = dialog('Name',title,...,'Position',FigPos);    
    'Visible'         ,'off'                      , ...
    'Name'            ,title                      , ...
    'Pointer'         ,'arrow'                    , ...
    'Position'        ,FigPos                     , ...
    'KeyPressFcn'     ,@doFigureKeyPress          , ...
    'IntegerHandle'   ,'off'                      , ...
    'WindowStyle'     ,'normal'                   , ...
    'HandleVisibility','callback'                 );%, ...
%     'CloseRequestFcn' ,@doDelete                  );%, ...
% %      );%'Tag'             ,Title                        ...


%make the buttons.
number_buttons = numel(options);

button_positions = [ones(number_buttons,1)*(box_width-button_size(1))/2,...
        ((number_buttons)*button_step:-button_step:number_buttons)',...
        repmat(button_size,number_buttons,1)];
                     
    
for x = 1:length(options)
    button(x) = uicontrol('Parent',d, 'Style', 'pushbutton',...
        'String', options(x), ...
        'Position', button_positions(x,:),...
        'Callback',@button_callback,...
        'KeyPressFcn',@pushbutton1_KeyPressFcn);
end
    uicontrol(button(1))
                     
%wrtie message to figure. 
%this is mostly copied from the help for wraptext.m

% Create it in Units of Pixels, 100 wide, 10 high
pos1 = [150 100 colw 50];
ht1 = uicontrol('Style','Text','Position',pos1,'Parent',d);
string1 = {message};
outstring1 = textwrap(ht1,string1);
% Reset Units of ht1 to Characters to use the result
set(ht1,'Units','characters')
mespos1 = get(ht1,'Position');
                     
% Set height the length of the outstring1 cell array + 1.
mespos1(4) = length(outstring1)+1;
    
%get the other positions in characters so all the spacing can be done
%relative to the size of the message string
set(button,'Units','characters');
butpos = get(button(1),'Position');
set(d,'Units','characters');
gif_pos = get(d,'Position');

%set the message box size
mespos1(1) = (butpos(1) + butpos(3)/2) - mespos1(3)/2;
mespos1(2) = butpos(2) + mespos1(4)+.75;
set(ht1,'String',outstring1,'Position',mespos1,'Units','characters')
    
%set the diagloge box height
gif_pos(4) = mespos1(2) +  mespos1(4) +1.25;
set(d,'Position',gif_pos)
    

%position the figure nively on the screen.
set(d,'Units','pixels')
% FigPos    = ;
set(d,'Position',getnicedialoglocation(get(d,'Position'), get(d,'Units')))
% FigPos    = get(d,'Position');
% FigPos    = getnicedialoglocation(FigPos, get(d,'Units'));
% set(d,'Position',FigPos)
set(d,'Visible','on');

% Wait for d to close before running to completion
uiwait(d);

    function button_callback(source,~)
        %           choice = source.String;
        % This code uses dot notation to get properties.
        % Dot notation runs in R2014b and later.
        % For R2014a and earlier:
        choice = get(source,'String');
        choice = choice{:};
        delete(gcf);
    end

    function pushbutton1_KeyPressFcn(hObject, ~, ~)
        button_callback(hObject);
    end
end


function figure_size = getnicedialoglocation(figure_size, figure_units)
% adjust the specified figure position to fig nicely over GCBF
% or into the upper 3rd of the screen

%  Copyright 1999-2010 The MathWorks, Inc.
%  $Revision: 1.1.6.5 $

% parentHandle = gcbf;
% convertData.destinationUnits = figure_units;
% if ~isempty(parentHandle)
%     % If there is a parent figure
%     convertData.hFig = parentHandle;
%     convertData.size = get(parentHandle,'Position');
%     convertData.sourceUnits = get(parentHandle,'Units');  
%     c = []; 
% else
%     % If there is no parent figure, use the root's data
%     % and create a invisible figure as parent
%     convertData.hFig = figure('visible','off');
%     convertData.size = get(0,'ScreenSize');
%     convertData.sourceUnits = get(0,'Units');
%     c = onCleanup(@() close(convertData.hFig));
% end

% Get the size of the dialog parent in the dialog units
container_size = get(0,'ScreenSize');
% hgconvertunits(convertData.hFig, convertData.size ,...
%     convertData.sourceUnits, convertData.destinationUnits, get(convertData.hFig,'Parent'));

% delete(c);

figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));


end