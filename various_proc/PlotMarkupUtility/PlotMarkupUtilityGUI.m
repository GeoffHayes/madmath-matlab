function varargout = PlotMarkupUtilityGUI(varargin)
% PLOTMARKUPUTILITYGUI MATLAB code for PlotMarkupUtilityGUI.fig
%      PLOTMARKUPUTILITYGUI, by itself, creates a new PLOTMARKUPUTILITYGUI or raises the existing
%      singleton*.
%
%      H = PLOTMARKUPUTILITYGUI returns the handle to a new PLOTMARKUPUTILITYGUI or the handle to
%      the existing singleton*.
%
%      PLOTMARKUPUTILITYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTMARKUPUTILITYGUI.M with the given input arguments.
%
%      PLOTMARKUPUTILITYGUI('Property','Value',...) creates a new PLOTMARKUPUTILITYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotMarkupUtilityGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotMarkupUtilityGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotMarkupUtilityGUI

% Last Modified by GUIDE v2.5 04-Jun-2014 20:04:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlotMarkupUtilityGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PlotMarkupUtilityGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before PlotMarkupUtilityGUI is made visible.
function PlotMarkupUtilityGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotMarkupUtilityGUI (see VARARGIN)

% Choose default command line output for PlotMarkupUtilityGUI
handles.output = hObject;

handles.markupUtil = PlotMarkupUtility([1 0 0 ]);

set(handles.txtMarkupColour,'BackgroundColor',[1 0 0]);

set(handles.sldrMarkupWidth,'SliderStep',[1 0.5],'Max',50,'Min',0.5,'Value',0.5);

if ~isfield(handles,'hListener')
    handles.hListener = ...
        addlistener(handles.sldrMarkupWidth,'ContinuousValueChange', ...
        @respondToContSlideCallback);
end

set(handles.txtMarkupWidth,'String',num2str(0.5));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlotMarkupUtilityGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlotMarkupUtilityGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnChangeMarkupColour.
function btnChangeMarkupColour_Callback(hObject, eventdata, handles)
% hObject    handle to btnChangeMarkupColour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

colour = get(handles.txtMarkupColour,'BackgroundColor');
colour = uisetcolor(colour);
set(handles.txtMarkupColour,'BackgroundColor',colour);
util = handles.markupUtil;
handles.markupUtil = util.setMarkupColour(colour);
guidata(hObject,handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles)
    if isfield(handles,'markupUtil')
        delete(handles.markupUtil);
    end
end

% --- Executes on slider movement.
function sldrMarkupWidth_Callback(hObject, eventdata, handles)
% hObject    handle to sldrMarkupWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sldrMarkupWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sldrMarkupWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function respondToContSlideCallback(hObject, eventdata)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% first we need the handles structure which we can get from hObject
handles = guidata(hObject);

if ~isempty(handles)
    if isfield(handles,'markupUtil')
        util = handles.markupUtil;
        val  = get(hObject,'Value');
        handles.markupUtil = util.setMarkupLineWidth(val);
        set(handles.txtMarkupWidth,'String',num2str(val));
        guidata(hObject,handles);
    end
end

% --------------------------------------------------------------------
function uipushtoolNewFig_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolNewFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles)
    if isfield(handles,'markupUtil')
        util = handles.markupUtil;
        handles.markupUtil = util.newFigure;
        guidata(hObject,handles);
    end
end

% --------------------------------------------------------------------
function uipushtoolUndoMarkup_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolUndoMarkup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles)
    if isfield(handles,'markupUtil')
        util = handles.markupUtil;
        handles.markupUtil = util.removeLastMarkup;
        guidata(hObject,handles);
    end
end
