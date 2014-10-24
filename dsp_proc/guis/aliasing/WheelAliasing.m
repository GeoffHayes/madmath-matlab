function varargout = WheelAliasing(varargin)
% WHEELALIASING MATLAB code for WheelAliasing.fig
%      WHEELALIASING, by itself, creates a new WHEELALIASING or raises the existing
%      singleton*.
%
%      H = WHEELALIASING returns the handle to a new WHEELALIASING or the handle to
%      the existing singleton*.
%
%      WHEELALIASING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WHEELALIASING.M with the given input arguments.
%
%      WHEELALIASING('Property','Value',...) creates a new WHEELALIASING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WheelAliasing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WheelAliasing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WheelAliasing

% Last Modified by GUIDE v2.5 24-Aug-2014 22:47:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WheelAliasing_OpeningFcn, ...
                   'gui_OutputFcn',  @WheelAliasing_OutputFcn, ...
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


% --- Executes just before WheelAliasing is made visible.
function WheelAliasing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WheelAliasing (see VARARGIN)

% Choose default command line output for WheelAliasing
handles.output = hObject;

% draw a four-spoked wheel on the axes
radius = 5;
hold(handles.axes1,'on');
handles.spoke1 = plot(handles.axes1,0:0.01:radius,zeros(1,length(0:0.01:radius)),'r','LineWidth',2);
handles.spoke2 = plot(handles.axes1,zeros(1,length(0:0.01:radius)),0:-0.01:-radius,'b','LineWidth',2);
handles.spoke3 = plot(handles.axes1,0:-0.01:-radius,zeros(1,length(0:0.01:radius)),'g','LineWidth',2);
handles.spoke4 = plot(handles.axes1,zeros(1,length(0:0.01:radius)),0:0.01:radius,'m','LineWidth',2);
xlim(handles.axes1,[-8,8]);
ylim(handles.axes1,[-8,8]);

% create the wheel
[C] = func_circle(0,0,radius);
handles.wheeltop = plot(handles.axes1,C(:,1),C(:,2),'k-');
handles.wheelbot = plot(handles.axes1,C(:,1),C(:,3),'k-');

handles.freq = 5;
handles.fs   = 5;

handles.timer  = [];

set(handles.frequency, 'String',num2str(handles.freq));
set(handles.sampleRate,'String',num2str(handles.fs));

% Update handles structure
guidata(hObject, handles);

set(handles.stop,'Enable','off');

% UIWAIT makes WheelAliasing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WheelAliasing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.timer)
    try
        stop(handles.timer);
        delete(handles.timer);
        handles.timer = [];
    catch
        % intentionally left blank
    end
end

% check the samples per second
fs = 1/str2num(char(get(handles.sampleRate,'String'))); 

handles.timer = timer('Name','WheelAliasing',                   ...
                      'Period',fs,                     ... 
                      'StartDelay',0,                 ... 
                      'TasksToExecute',inf,           ... 
                      'ExecutionMode','fixedSpacing', ...
                      'TimerFcn',{@timerCallback,handles.figure1}); 
guidata(hObject,handles);

start(handles.timer);

set(handles.sampleRate,'Enable','off');
set(handles.play,'Enable','off');
set(handles.stop,'Enable','on');

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    stop(handles.timer);
    delete(handles.timer);
    handles.timer = [];
    guidata(hObject,handles);
catch
    % intentionally left blank
end

set(handles.sampleRate,'Enable','on');
set(handles.play,'Enable','on');
set(handles.stop,'Enable','off');



function frequency_Callback(hObject, eventdata, handles)
% hObject    handle to frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency as text
%        str2double(get(hObject,'String')) returns contents of frequency as a double

val = str2num(char(get(hObject,'String')));

if ~isempty(val)
    if val>500
        set(hObject,'String','500');
    end
else
    set(hObject,'String','5');
end

handles.freq = str2num(char(get(hObject,'String')));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to sampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampleRate as text
%        str2double(get(hObject,'String')) returns contents of sampleRate as a doubleval = str2num(char(get(hObject,'String')));
val = str2num(char(get(hObject,'String')));

if ~isempty(val)
    if val>75
        set(hObject,'String','75');
    end
else
    set(hObject,'String','5');
end

handles.fs = str2num(char(get(hObject,'String')));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [] = timerCallback(~,~,guiHandle)

try
    if ~isempty(guiHandle)        
        handles = guidata(guiHandle);

        numRevols = handles.freq*1/handles.fs;
        
        rDegs = 360.0*numRevols;
        
        R = [cosd(rDegs) sind(rDegs); -sind(rDegs) cosd(rDegs)];
        
        % rotate the spokes
        xdata1 = get(handles.spoke1,'XData');
        ydata1 = get(handles.spoke1,'YData');
        
        xdata2 = get(handles.spoke2,'XData');
        ydata2 = get(handles.spoke2,'YData');
        
        xdata3 = get(handles.spoke3,'XData');
        ydata3 = get(handles.spoke3,'YData');
        
        xdata4 = get(handles.spoke4,'XData');
        ydata4 = get(handles.spoke4,'YData');
        
        for k=1:length(xdata1)
            u = R*[xdata1(k);ydata1(k)];
            xdata1(k) = u(1);
            ydata1(k) = u(2);
            
            u = R*[xdata2(k);ydata2(k)];
            xdata2(k) = u(1);
            ydata2(k) = u(2);
            
            u = R*[xdata3(k);ydata3(k)];
            xdata3(k) = u(1);
            ydata3(k) = u(2);
            
            u = R*[xdata4(k);ydata4(k)];
            xdata4(k) = u(1);
            ydata4(k) = u(2);
        end
        
        set(handles.spoke1,'XData',xdata1,'YData',ydata1);
        set(handles.spoke2,'XData',xdata2,'YData',ydata2);
        set(handles.spoke3,'XData',xdata3,'YData',ydata3);
        set(handles.spoke4,'XData',xdata4,'YData',ydata4);
        

    end
catch me
    % intentionally left blank
    me
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    stop(handles.timer);
    delete(handles.timer);
    handles.timer = [];
    guidata(hObject,handles);
catch
    % intentionally left blank
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
selection = questdlg('Close Wheel Aliasing?',...
                     'Close Request Function',...
                     'Yes','No','Yes');
switch selection,
   case 'Yes',
     delete(hObject);
   case 'No'
     return
end
