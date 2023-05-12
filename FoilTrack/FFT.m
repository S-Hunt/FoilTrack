function varargout = FFT(varargin)
% FFT MATLAB code for FFT.fig
%      FFT, by itself, creates a new FFT or raises the existing
%      singleton*.
%
%      H = FFT returns the handle to a new FFT or the handle to
%      the existing singleton*.
%
%      FFT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FFT.M with the given input arguments.
%
%      FFT('Property','Value',...) creates a new FFT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FFT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FFT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FFT

% Last Modified by GUIDE v2.5 18-Oct-2019 13:20:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FFT_OpeningFcn, ...
                   'gui_OutputFcn',  @FFT_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
    %imshow(axes1, varargin{1})
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before FFT is made visible.
function FFT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FFT (see VARARGIN)

% Choose default command line output for FFT
handles.output = hObject

%get the image to porcess
if isstr(varargin{1})
    handles.image = imread(varargin{1});
    handles.image_name = varargin{1};
elseif nargin >=2
    handles.image = varargin{1};
    handles.image_name = varargin{2};
else
    error('Input is not recognised')
end
handles.FilterIm = logical(zeros(size(handles.image)));

% Update handles structure
guidata(hObject, handles)

% This sets up the initial plot - only do when we are invisible
% so window can get raised using FFT.
if strcmp(get(hObject,'Visible'),'off')
    
    imagesc(handles.axes1, handles.image)
    
    Fimage = fft2(handles.image);
    Fplot  = log10(fftshift(abs(Fimage)));

    imagesc(handles.axes2, Fplot)
    colormap(handles.axes2, 'gray')
    handles.axes2.CLim = [min(Fplot(:)), max(Fplot(:))];
    
    handles.editMAX.String = sprintf('%3.2g',max(Fplot(:)));
    handles.editMin.String = sprintf('%3.2g',min(Fplot(:)));
end


handles.filter_state = 0;

%link the axes in the image plots.
linkaxes([handles.axes1 handles.axes3 handles.axes4],'xy')


process_Callback(hObject, eventdata, handles)



% UIWAIT makes FFT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FFT_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)
% hObject    handle to process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Fim = fft2(handles.image);
 %Fplot  = log10((fftshift(abs(Fimage))));

 %make image mask from handle.filter.
 filt = zeros(size(handles.image));
 if isfield(handles, 'filter')
     for x = 1:length(handles.filter)
        filter_t = createMask(handles.filter{x}, handles.axes3.Children);
        filt = filt+filter_t;
     end
     
 end
 if handles.filter_state == 0
     filt = ~filt;
 end
%  if handles.filter_state == 0
%     Fimage_filtered = fftshift(Fim).*~filt;
%     handles.FilterIm = ~filt;
%  else
%     Fimage_filtered = fftshift(Fim).*filt;
%     handles.FilterIm = filt;
%  end
%  
%  filtered_im = real(ifft2(ifftshift(Fimage_filtered)));

 filtered_im = FFTImageFilter(handles.image, filt);
 handles.FilterIm = filt;
 
 imagesc(handles.axes3, filtered_im)
 imagesc(handles.axes4, double(handles.image)-filtered_im)
 colormap(handles.axes4, 'cool')
 
 % Update handles structure
guidata(hObject, handles)
 


% --- Executes on button press in pushbuttonPoly.
function pushbuttonPoly_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = impoly(handles.axes2);

if isfield(handles, 'filter')
    handles.filter{end+1} = h
else
    handles.filter{1} = h
end

%keyboard
% Update handles structure
guidata(hObject, handles)



% --- Executes on button press in pushbuttonEllipse.
function pushbuttonEllipse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEllipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = imellipse(handles.axes2)

if isfield(handles, 'filter')
    handles.filter{end+1} = h
else
    handles.filter{1} = h
end

% Update handles structure
guidata(hObject, handles)


% --- Executes on button press in pushbuttonRect.
function pushbuttonRect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = imrect(handles.axes2)

if isfield(handles, 'filter')
    handles.filter{end+1} = h
else
    handles.filter{1} = h
end

% Update handles structure
guidata(hObject, handles)

% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%pushbuttonSave_Callback(hObject, eventdata, handles)

delete(handles.figure1)



% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get 'base' file name
[~,fname,~] = fileparts(handles.image_name);

if isfield(handles, 'filter')
    flt = handles.filter;
    flt_state = handles.filter_state;
    
    save_vars = {'flt','flt_state'};
    whos flt flt_state
    save([fname,'.mat'], save_vars{:})
    
end

figure
imagesc(handles.FilterIm)
colormap gray

imwrite(handles.FilterIm,[fname,'_filter.png'])



% --- Executes on button press in pushbuttonRotate.
function pushbuttonRotate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

im_size = size(handles.image);

for x = 1: length(handles.filter)
    Pos = handles.filter{x}.getPosition;
    

    if isa(handles.filter{x}, 'impoly')
        NewPos = [im_size(2)-Pos(:,1), im_size(1)-Pos(:,2)];
        h = impoly(handles.axes2, NewPos);
    elseif isa(handles.filter{x}, 'imellipse')
        bot = Pos(2);
        top = Pos(2) + Pos(4);
        new_bot = im_size(1) - top;
        l = Pos(1);
        r = Pos(1) + Pos(3);
        new_l = im_size(2) - r;
        %new_top = im_size(1) - bot;
        NewPos = [new_l, new_bot Pos(3), Pos(4)];
        h = imellipse(handles.axes2, NewPos);
    elseif isa(handles.filter{x}, 'imrect')
        bot = Pos(2);
        top = Pos(2) + Pos(4);
        new_bot = im_size(1) - top;
        l = Pos(1);
        r = Pos(1) + Pos(3);
        new_l = im_size(2) - r;
        %new_top = im_size(1) - bot;
        NewPos = [new_l, new_bot Pos(3), Pos(4)];
        h = imrect(handles.axes2, NewPos);
    else
        error 'oops'
    end
    handles.filter{end+1} = h;
end

% Update handles structure
guidata(hObject, handles)


% --- Executes on button press in pushbuttonMirrorTB.
function pushbuttonMirrorTB_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMirrorTB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

im_size = size(handles.image);

for x = 1: length(handles.filter)
    Pos = handles.filter{x}.getPosition;
    

    if isa(handles.filter{x}, 'impoly')
        NewPos = [Pos(:,1), im_size(1)-Pos(:,2)];
        h = impoly(handles.axes2, NewPos);
    elseif isa(handles.filter{x}, 'imellipse')
        bot = Pos(2);
        top = Pos(2) + Pos(4);
        new_bot = im_size(1) - top;
        %new_top = im_size(1) - bot;
        NewPos = [Pos(1), new_bot, Pos(3), Pos(4)];
        h = imellipse(handles.axes2, NewPos);
    elseif isa(handles.filter{x}, 'imrect')
        bot = Pos(2);
        top = Pos(2) + Pos(4);
        new_bot = im_size(1) - top;
        %new_top = im_size(1) - bot;
        NewPos = [Pos(1), new_bot, Pos(3), Pos(4)];
        h = imrect(handles.axes2, NewPos);
    else
        error 'oops'
    end
    handles.filter{end+1} = h;
end

%keyboard
% Update handles structure
guidata(hObject, handles)



% --- Executes on button press in pushbuttonMirrorLR.
function pushbuttonMirrorLR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMirrorLR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

im_size = size(handles.image);

for x = 1: length(handles.filter)
    Pos = handles.filter{x}.getPosition;
    

    if isa(handles.filter{x}, 'impoly')
        NewPos = [im_size(2)-Pos(:,1), Pos(:,2)];
        h = impoly(handles.axes2, NewPos);
    elseif isa(handles.filter{x}, 'imellipse')
        l = Pos(1);
        r = Pos(1) + Pos(3);
        new_bot = im_size(2) - r;
        %new_top = im_size(1) - bot;
        NewPos = [new_bot, Pos(2) Pos(3), Pos(4)];
        h = imellipse(handles.axes2, NewPos);
    elseif isa(handles.filter{x}, 'imrect')
        l = Pos(1);
        r = Pos(1) + Pos(3);
        new_bot = im_size(2) - r;
        %new_top = im_size(1) - bot;
        NewPos = [new_bot, Pos(2) Pos(3), Pos(4)];
        h = imrect(handles.axes2, NewPos);
    else
        error 'oops'
    end
    handles.filter{end+1} = h;
end

%keyboard
% Update handles structure
guidata(hObject, handles)



% --- Executes on button press in checkbox1. %filter inside or outside.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

if get(hObject,'Value')
    handles.filter_state = 1;
else
    handles.filter_state = 0;
end

% Update handles structure
guidata(hObject, handles)
    



function editMAX_Callback(hObject, eventdata, handles)
% hObject    handle to editMAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMAX as text
%        str2double(get(hObject,'String')) returns contents of editMAX as a double

lims = get(handles.axes2, 'Clim');
set(handles.axes2, 'Clim', [lims(1), str2double(get(hObject,'String'))]);



% --- Executes during object creation, after setting all properties.
function editMAX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMin_Callback(hObject, eventdata, handles)
% hObject    handle to editMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMin as text
%        str2double(get(hObject,'String')) returns contents of editMin as a double

lims = get(handles.axes2, 'Clim');
set(handles.axes2, 'Clim', [str2double(get(hObject,'String')), lims(2)]);


% --- Executes during object creation, after setting all properties.
function editMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when uipanel1 is resized.
function uipanel1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
