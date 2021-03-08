function Answer = SelectBoxes_tableUI(varargin)
% function to make table for positioning of individual boxes.
% Part of SelectBoxes
% Simon Hunt August 2017

%parse inputs
if nargin  < 2
    error('Not enough inputs')
end
fit_opts = varargin{1};%SelectBoxes_position('posibilities'); %call possibilities from SelectBoxes_position
BoxX = varargin{2};%[ 1 2; 3 4; 5 6 ; 7 8 ; 9 0; 1 2 ];
Title = 'Box position setting';

%make table text
if size(BoxX,1)<=10
    dat =  [num2cell((1:size(BoxX,1))') repmat({' '},size(BoxX,1),1)];%   {1, ''; 2, ''; 3, ''; 4, ''; 5, ''; 6, ''; 7, ''; 8, ''; 9, ''; 10,''};
else
    error('SelectBoxes_tableUI cannot deal with more than 10 boxes. Fit by row option is not created');
end
columnname =   {'Box Number', 'Fitting method'};
columnformat = {'numeric', fit_opts'};
columneditable =  [true true];

% Setup button height and width
BtnWidth  = 53;%max(DefBtnWidth,BtnExtent(3)+8);
BtnHeight = 23;%max(DefBtnHeight,BtnExtent(4)*btnMargin);
DefOffset    = 5;

FigWidth=175;
FigHeight=100;
FigPos(3:4)=[FigWidth FigHeight];  %#ok
FigColor=get(0,'DefaultUicontrolBackgroundColor');

%make the figure and table ui
h=dialog(...
    'KeyPressFcn'      ,@doFigureKeyPress, ...
    'Name'             ,Title      , ...
    'Pointer'          ,'arrow'    , ...
    'Units'            ,'pixels'   , ...
    'UserData'         ,'Cancel'   , ...
    'Tag'              ,Title      , ...
    'HandleVisibility' ,'callback' , ...
    'Color'            ,FigColor   );%, ...;

t = uitable(h, 'Units','normalized','Position',...
    [0.1 0.1 0.8 0.8], 'Data', dat,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', columneditable,...
    'ColumnWidth',{100, 250},...
    'RowName',[]);
set(h,'CloseRequestFcn',@myCloseFcn)
%give a unique Tag:
set(h,'Tag', 'myTag')
set(t,'Tag','myTableTag')


% add buttons
OKHandle=uicontrol(h     ,              ...  BtnInfo      , ...
  'Position'   ,[ FigWidth-2*BtnWidth-2*DefOffset DefOffset BtnWidth BtnHeight ] , ...
  'KeyPressFcn',@doControlKeyPress , ...
  'String'     ,getString(message('MATLAB:uistring:popupdialogs:OK'))        , ...
  'Callback'   ,@doCallback , ...
  'Tag'        ,'OK'        , ...
  'UserData'   ,'OK'          ...
  );
CancelHandle=uicontrol(h     ,              ...  BtnInfo      , ...
  'Position'   ,[ FigWidth-BtnWidth-DefOffset DefOffset BtnWidth BtnHeight ]           , ...
  'KeyPressFcn',@doControlKeyPress            , ...
  'String'     ,getString(message('MATLAB:uistring:popupdialogs:Cancel'))    , ...
  'Callback'   ,@doCallback , ...
  'Tag'        ,'Cancel'    , ...
  'UserData'   ,'Cancel'       ...
  ); %#ok


uiwait(h);

% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
if ishghandle(h)
  if strcmp(get(h,'UserData'),'OK'),
    Answer=get(t, 'Data');
    
    %FIXME: This will only work for individual boxes and not for rows.
    Answer = Answer(:,2);
  end
  delete(h);
else
  Answer={};
end


end

function myCloseFcn(~,~)
myfigure=findobj('Tag','myTag');
myData=get(findobj(myfigure,'Tag','myTableTag'),'Data');
assignin('base','myTestData',myData)
delete(myfigure)
end

function doControlKeyPress(obj, evd) %#ok
switch(evd.Key)
  case {'return'}
    if ~strcmp(get(obj,'UserData'),'Cancel')
      set(gcbf,'UserData','OK');
      uiresume(gcbf);
    else
      delete(gcbf)
    end
  case 'escape'
    delete(gcbf)
end
% myCloseFcn
end

function doCallback(obj, evd) %#ok
if ~strcmp(get(obj,'UserData'),'Cancel')
  set(gcbf,'UserData','OK');
  uiresume(gcbf);
else
  delete(gcbf)
end
% myCloseFcn
end