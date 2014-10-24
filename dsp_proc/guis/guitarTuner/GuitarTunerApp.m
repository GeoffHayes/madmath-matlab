function varargout = GuitarTunerApp(varargin)
% GUITARTUNERAPP MATLAB code for GuitarTunerApp.fig
%      GUITARTUNERAPP, by itself, creates a new GUITARTUNERAPP or raises the existing
%      singleton*.
%
%      H = GUITARTUNERAPP returns the handle to a new GUITARTUNERAPP or the handle to
%      the existing singleton*.
%
%      GUITARTUNERAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITARTUNERAPP.M with the given input arguments.
%
%      GUITARTUNERAPP('Property','Value',...) creates a new GUITARTUNERAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuitarTunerApp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuitarTunerApp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuitarTunerApp

% Last Modified by GUIDE v2.5 29-Sep-2014 11:39:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuitarTunerApp_OpeningFcn, ...
                   'gui_OutputFcn',  @GuitarTunerApp_OutputFcn, ...
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


% --- Executes just before GuitarTunerApp is made visible.
function GuitarTunerApp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuitarTunerApp (see VARARGIN)

% Choose default command line output for GuitarTunerApp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

resetAppData(hObject);

% create the audio recorder
handles = guidata(hObject);
handles.recorder            = audiorecorder(handles.Fs,8,1);
set(handles.recorder,'TimerPeriod',1,'TimerFcn',{@audioTimer,hObject});

guidata(hObject, handles);

set(handles.uitoggletoolPlay,'Enable','on');
set(handles.uitoggletoolPause,'Enable','off');
set(handles.uitoggletoolStop,'Enable','off');


% UIWAIT makes GuitarTunerApp wait for user response (see UIRESUME)
% uiwait(handles.figGuitarTuner);


% --- Outputs from this function are returned to the command line.
function varargout = GuitarTunerApp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function uitoggletoolPlay_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.uitoggletoolPlay,'Enable','off');
set(handles.uitoggletoolPause,'Enable','on');
set(handles.uitoggletoolStop,'Enable','on');

resetAppData(handles.figGuitarTuner);

record(handles.recorder);


% --------------------------------------------------------------------
function uitoggletoolPause_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isrecording(handles.recorder)
    pause(handles.recorder);
else
    resume(handles.recorder);
end


% --------------------------------------------------------------------
function uitoggletoolStop_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stop(handles.recorder);

set(handles.uitoggletoolPlay,'Enable','on','State','off');
set(handles.uitoggletoolPause,'Enable','off','State','off');
set(handles.uitoggletoolStop,'Enable','off','State','off');

function resetAppData(hObject)

% get the handles structure
handles = guidata(hObject);

handles.Fs                  = 8192;
handles.N                   = 8192;

handles.dispFreqMaxHz       = 800;

handles.lastSampleAt        = 0;

handles.spectrogramImgQntzd = uint8(zeros(handles.N,512,3));
handles.spectrogramImgDb    = zeros(handles.N,512);
handles.spectrogramImgCol   = 1;

handles.minDb               = -90;
handles.maxDb               = 100;
handles.dbLevels            = linspace(handles.minDb,handles.maxDb,255);
handles.spectrogramClrMap   = uint8(255*colormap(hot(255)));

set(handles.radiobuttonLowE,'Value',1);

handles.stringFrequencies = [82.41 110.00 146.83 196.00 246.94 329.63];
handles.targStrFreq       = handles.stringFrequencies(1);

% save the data
guidata(hObject,handles);

function audioTimer(hObject,varargin)

hFigure = varargin{2};
handles = guidata(hFigure);

samples  = getaudiodata(hObject);

Fs       = handles.Fs;

% skip if not enough data
if length(samples)<handles.lastSampleAt+1+Fs
    return;
end

X   = samples(handles.lastSampleAt+1:handles.lastSampleAt+1+Fs);
Y   = r2fft(X,Fs);
dbY = 20*log10(abs(Y));
dbQ = repmat(handles.spectrogramClrMap(end,:),length(dbY),1);

for k=1:length(dbY)
    for u=25:length(handles.dbLevels)
        if handles.dbLevels(u)>dbY(k)
            dbQ(k,:) = handles.spectrogramClrMap(u,:);
            break;
        end
    end
end

if handles.spectrogramImgCol>size(handles.spectrogramImgQntzd,2)
    handles.spectrogramImgQntzd(:,1:end-1,:) = handles.spectrogramImgQntzd(:,2:end,:);
    handles.spectrogramImgQntzd(:,end,:)     = dbQ;
    handles.spectrogramImgDb(:,1:end-1)      = handles.spectrogramImgDb(:,2:end);
    handles.spectrogramImgDb(:,end)          = dbY;
    handles.timeOffsetSecs                   = handles.spectrogramImgCol - size(handles.spectrogramImgQntzd,2) + 1;
else
    handles.spectrogramImgQntzd(:,handles.spectrogramImgCol,:) = dbQ;
    handles.spectrogramImgDb(:,handles.spectrogramImgCol)      = dbY;
    handles.timeOffsetSecs                                     = 0;
end

image(handles.spectrogramImgQntzd(1:handles.dispFreqMaxHz,:,:),'Parent',handles.axesSpectrogram);
handles.spectrogramImgCol = handles.spectrogramImgCol + 1;

set(handles.axesSpectrogram,...
    'YTick',linspace(1,handles.dispFreqMaxHz,8),...
    'YTickLabel',num2str(floor(linspace(1,handles.dispFreqMaxHz,8)')), ...        
    'XTick',linspace(1,size(handles.spectrogramImgQntzd,2),5),...
    'XTickLabel',num2str(floor(linspace(1,size(handles.spectrogramImgQntzd,2),5))'+handles.timeOffsetSecs), ...                      
    'Tag','axesSpectrogram');

handles.lastSampleAt = handles.lastSampleAt + Fs;

guidata(hFigure,handles);


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
