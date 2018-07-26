function varargout = threshPick(varargin)
%by: Joey Broussard
%at UC Davis 09/08/14
%
%Used in conjunction with maskUpdate which acts as a function callback for the
%GUI.
%This function calls a GUI which displays allows a user to adjust the level
%of BW masking. Input is an image for thresholding. Output is the
%thresholded image mask.


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @threshPick_OpeningFcn, ...
                   'gui_OutputFcn',  @threshPick_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before threshPick is made visible.
function threshPick_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output
handles.output = hObject;

%
global Img
Img = varargin{1};
assignin('base','myGUIvar',Img);
axes(handles.axes1);
imshow(Img);


handles.output = Img;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes threshPick wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = threshPick_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see varargout);

global mask
global level

% Get default command line output from handles structure
varargout{1} = mask;
varargout{2} = level;


% --- Executes on slider movement.
function Level_Callback(hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

maskUpdate(handles);

% --- Executes during object creation, after setting all properties.
function Level_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in EndButton.
function EndButton_Callback(hObject, eventdata, handles)

close(handles.figure1)


function Display_Callback(hObject, eventdata, handles)

maskUpdate(handles);

% --- Executes during object creation, after setting all properties.
function Display_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
