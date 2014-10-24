function varargout = Radix2FFTApp(varargin)
% RADIX2FFTAPP MATLAB code for Radix2FFTApp.fig
%      RADIX2FFTAPP, by itself, creates a new RADIX2FFTAPP or raises the existing
%      singleton*.
%
%      H = RADIX2FFTAPP returns the handle to a new RADIX2FFTAPP or the handle to
%      the existing singleton*.
%
%      RADIX2FFTAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RADIX2FFTAPP.M with the given input arguments.
%
%      RADIX2FFTAPP('Property','Value',...) creates a new RADIX2FFTAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Radix2FFTApp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Radix2FFTApp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Radix2FFTApp

% Last Modified by GUIDE v2.5 17-Sep-2014 23:45:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Radix2FFTApp_OpeningFcn, ...
                   'gui_OutputFcn',  @Radix2FFTApp_OutputFcn, ...
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


% --- Executes just before Radix2FFTApp is made visible.
function Radix2FFTApp_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Radix2FFTApp (see VARARGIN)

% Choose default command line output for Radix2FFTApp
handles.output = hObject;

handles.timer  = [];

% create the listener for the max dB slider
handles.sliderMaxDbListener = ...
    addlistener(handles.sliderMaxDb,'ContinuousValueChange', ...
                @(guiHandle,eventdata) sliderMaxDbContValCallback(...
                hObject,eventdata));
            
handles.sliderSpecMaxFreqListener = ...
    addlistener(handles.sliderSpecMaxFreqHz,'ContinuousValueChange', ...
                @(guiHandle,eventdata) sliderSpecMaxFreqContValCallback(...
                hObject,eventdata));

% Update handles structure
guidata(hObject, handles);

% set the app data
setappdata(hObject,'sampleRateHz',512);
setappdata(hObject,'blockSize',256);
setappdata(hObject,'sineAmplitude',1);
setappdata(hObject,'sineFreqHz',5);
setappdata(hObject,'signalAmps',[]);
setappdata(hObject,'signalFreqsHz',[]);
setappdata(hObject,'addNoise',false);
setappdata(hObject,'displayOn',true(6,1));
setappdata(hObject,'inputfcnHandle',[]);

set(handles.checkboxAnalogData,'Value',1);
set(handles.checkboxDigitalData,'Value',1);
set(handles.checkboxFFTData,'Value',1);
set(handles.checkboxPhaseData,'Value',1);
set(handles.checkboxRealData,'Value',1);
set(handles.checkboxImagData,'Value',1);

set(handles.uitoggletoolStart,'Enable','off');

resetAppData(hObject);


% UIWAIT makes Radix2FFTApp wait for user response (see UIRESUME)
% uiwait(handles.figRadix2FFT);


% --- Outputs from this function are returned to the command line.
function varargout = Radix2FFTApp_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function uipushtoolBlockSize_ClickedCallback(hObject, ~, handles)
% hObject    handle to uipushtoolBlockSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt    = 'Radix-2 block size: ';
title     = 'N Input';
numLines  = 1;
hfig      = handles.figRadix2FFT;
currValue = {num2str(getappdata(hfig,'blockSize'))};
userInput = inputdlg(prompt,title,numLines,currValue);

if ~isempty(userInput)
    newBlockSize = str2double(char(userInput));
    if ~isempty(userInput)
        newBlockSize = 2^nextpow2(ceil(newBlockSize));
        setappdata(hfig,'blockSize',newBlockSize);
        setappdata(hfig,'spectrogramImgQntzd',uint8(zeros(newBlockSize,512,3)));
        setappdata(hObject,'spectrogramImgDb',zeros(newBlockSize,512));
        setappdata(hfig,'spectrogramImgCol',1);
    end
end


% --------------------------------------------------------------------
function uipushtoolSamplingRate_ClickedCallback(~, ~, handles)
% hObject    handle to uipushtoolSamplingRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt    = 'Sample Rate (Hz): ';
title     = 'Fs Input';
numLines  = 1;
hfig      = handles.figRadix2FFT;
currValue = {num2str(getappdata(hfig,'sampleRateHz'))};
userInput = inputdlg(prompt,title,numLines,currValue);

if ~isempty(userInput)
    newSampleRate = str2double(char(userInput));
    if ~isempty(userInput)
        newSampleRate = ceil(abs(newSampleRate));
        setappdata(hfig,'sampleRateHz',newSampleRate);
        set(handles.sliderSpecMaxFreqHz,'Max',newSampleRate,'Value',...
            min(newSampleRate,get(handles.sliderSpecMaxFreqHz,'Value')));
    end
end


% --------------------------------------------------------------------
function uipushtoolSignal_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if 1==0
    prompt    = {'Amplitude: ','Frequency (Hz): '};
    title     = 'Signal Input';
    numLines  = 1;
    hfig      = handles.figRadix2FFT;
    currValue = {num2str(getappdata(hfig,'sineAmplitude')),num2str(getappdata(hfig,'sineFreqHz'))};
    userInput = inputdlg(prompt,title,numLines,currValue);

    if ~isempty(userInput)
        newSineAmplitude = str2double(char(userInput{1}));
        newSineFreqHz    = str2double(char(userInput{2}));
        if ~isempty(userInput)
            setappdata(hfig,'sineAmplitude',abs(newSineAmplitude));
            setappdata(hfig,'sineFreqHz',abs(newSineFreqHz));
        end
    end
else
    [filename,~,filteridx] = uigetfile('*.m','Choose input function');

    if filteridx==1
        hfig      = handles.figRadix2FFT;
        setappdata(hfig,'inputfcnHandle',str2func(filename(1:end-2)));
        set(handles.uitoggletoolStart,'Enable','on');
    end
end

% --------------------------------------------------------------------
function uitoggletoolStart_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitoggletoolStart,'Enable','off');
set(handles.uitoggletoolStop,'Enable','on');

if ~isempty(handles.timer)
    try
        stop(handles.timer);
        delete(handles.timer);
        handles.timer = [];
    catch
        % intentionally left blank
    end
end

handles.timer = timer('Name','DSP',                   ...
                      'Period',0.93,                     ... 
                      'StartDelay',1,                 ... 
                      'TasksToExecute',inf,           ... 
                      'ExecutionMode','fixedSpacing', ...
                      'TimerFcn',{@timerCallback,handles.figRadix2FFT}); 
guidata(hObject,handles);

resetAppData(handles.figRadix2FFT);

% clear any persistent variables that may exist for the function
eval(['clear ' func2str(getappdata(handles.figRadix2FFT,'inputfcnHandle'))]);

start(handles.timer);

% --------------------------------------------------------------------
function uitoggletoolStop_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitoggletoolStart,'Enable','on','State','off');
set(handles.uitoggletoolStop,'Enable','off','State','off');

try
    stop(handles.timer);
    delete(handles.timer);
    handles.timer = [];
    guidata(hObject,handles);
catch
    % intentionally left blank
end

function [] = timerCallback(~,~,guiHandle)

try
    %h = findobj('Tag','figRadix2FFT');
    if ~isempty(guiHandle)
        
        handles = guihandles(guiHandle);
 
        t               = getappdata(guiHandle,'timestampsSecs');
        fs              = getappdata(guiHandle,'sampleRateHz');
        N               = getappdata(guiHandle,'blockSize');
        Dunused         = getappdata(guiHandle,'digTimeDomainData');
        addNoise        = getappdata(guiHandle,'addNoise');
        displayOn       = getappdata(guiHandle,'displayOn');
        fcnHandle       = getappdata(guiHandle,'inputfcnHandle');
        specImgQntzd    = getappdata(guiHandle,'spectrogramImgQntzd');
        specImgDb       = getappdata(guiHandle,'spectrogramImgDb');
        specImgCol      = getappdata(guiHandle,'spectrogramImgCol');
        specImgDbLevels = getappdata(guiHandle,'spectrogramDbLevels');
        specImgClrMap   = getappdata(guiHandle,'spectrogramClrMap');

        % create the new timestamp vector
        if isempty(t)
            maxTimestampSec = 0;
        else
            maxTimestampSec = max(t);
        end
        
        % use a gap of 1/1000th of a second
        analogGap = 0.00001;
        
        numSamplesPerSecond = 1/analogGap;
        if fs > numSamplesPerSecond
            fs = numSamplesPerSecond;
        end    
        ts = 1/fs;
        
        td  = (0:ts:1-ts)';
        tc  = (0:analogGap:td(end))' + maxTimestampSec;
        td  = td + maxTimestampSec;
        tdi = floor(linspace(1,length(tc),length(td)));
        
        maxTimestampSec = td(end);

%         amp  = getappdata(guiHandle,'sineAmplitude');
%         frq  = getappdata(guiHandle,'sineFreqHz');
 
        %A = amp*sin(2*pi*frq*tc) + noise;
        A = fcnHandle(tc);
        
        if length(A) ~= length(tc)
           tc  = linspace(tc(1),tc(end),length(A)); 
           tdi = floor(linspace(1,length(tc),length(td)));
        end
        
        if addNoise
            A = A + addNoise*randn(size(tc));
        end
        
        D = A(tdi);
        
        if displayOn(1)
            plot(handles.axesAnalog,tc,A,'g');
            line([tc(1) tc(end)],[0 0],'Color','w','Parent',handles.axesAnalog);
            set(handles.axesAnalog,'Color','k','Tag','axesAnalog');
            maxAbsA = max(abs(min(A)),abs(max(A)));
            axis(handles.axesAnalog,[tc(1) tc(end) -maxAbsA maxAbsA]);
        else
            cla(handles.axesAnalog);
            set(handles.axesAnalog,'Color','k');
        end

        if displayOn(2)
            stem(handles.axesDigital,td,D,'g','MarkerSize',2,'MarkerFaceColor','g');
            line([td(1) td(end)],[0 0],'Color','w','Parent',handles.axesDigital);
            set(handles.axesDigital,'Color','k','Tag','axesDigital');
            maxAbsD = max(abs(min(D)),abs(max(D)));
            axis(handles.axesDigital,[td(1) td(end) -maxAbsD maxAbsD]);
        else
            cla(handles.axesDigital);
            set(handles.axesDigital,'Color','k');
        end
        
        % prepend leftover data from prior run
        if length(t)>1
            td = [t ; td];
            D = [Dunused ; D];
        end
        
        axesHandles = [handles.axesMagnitude,handles.axesReal,...
                       handles.axesPhase, handles.axesImag];
                   
        if ~displayOn(3)
            cla(handles.axesMagnitude);
            set(handles.axesMagnitude,'Color','k');
        end
        if ~displayOn(4)
            cla(handles.axesReal);
            set(handles.axesReal,'Color','k');
        end
        if ~displayOn(5)
            cla(handles.axesPhase);
            set(handles.axesPhase,'Color','k');
        end
        if ~displayOn(6)
            cla(handles.axesImag);
            set(handles.axesImag,'Color','k');
        end
                            
        axesHandles(~displayOn(3:end)) = -1;
        
        % do the radix-2 FFT
        doSpec = true;
        dispFreqMaxHz = floor(get(handles.sliderSpecMaxFreqHz,'Value'));
        while length(D)>=N
            Dp = D(1:N);
            Y  = r2fft(Dp,N);
            r2fftstats(Y,fs,Dp,axesHandles,dispFreqMaxHz);
            D  = D(N+1:end);
            td = td(N+1:end);
            
            if doSpec
                % convert to dB
                dbY = 20*log10(abs(Y));
                
                % quantize the data
                dbQ = uint8(zeros(length(dbY),1,3));
                
                for k=1:length(dbY)
                   idx = find(specImgDbLevels>dbY(k),1);
                   if isempty(idx)
                       dbQ(k,1,:) = specImgClrMap(end,:);
                   else
                       dbQ(k,1,:) = specImgClrMap(idx,:);
                   end
                end

                if specImgCol>size(specImgQntzd,2)
                    specImgQntzd(:,1:end-1,:) = specImgQntzd(:,2:end,:);
                    specImgQntzd(:,end,:)     = dbQ;
                    specImgDb(:,1:end-1)      = specImgDb(:,2:end);
                    specImgDb(:,end)          = dbY;
                    timeOffsetSecs            = specImgCol - size(specImgQntzd,2) + 1;
                else
                    specImgQntzd(:,specImgCol,:) = dbQ;
                    specImgDb(:,specImgCol)      = dbY;
                    timeOffsetSecs               = 0;
                end

                image(specImgQntzd(1:dispFreqMaxHz,:,:),'Parent',handles.axesSpectrogram);
                specImgCol = specImgCol + 1;
                setappdata(guiHandle,'spectrogramImgQntzd',specImgQntzd);
                setappdata(guiHandle,'spectrogramImgCol',specImgCol);
                setappdata(guiHandle,'spectrogramImgDb',specImgDb);

                set(handles.axesSpectrogram,...
                    'YTick',linspace(1,dispFreqMaxHz,8),...
                    'YTickLabel',num2str(floor(linspace(1,dispFreqMaxHz,8)')), ...        
                    'XTick',linspace(1,size(specImgQntzd,2),5),...
                    'XTickLabel',num2str(floor(linspace(1,size(specImgQntzd,2),5))'+timeOffsetSecs), ...                      
                    'Tag','axesSpectrogram');
                
                doSpec = false;
            end
            
        end
        
        set(handles.axesMagnitude,'Color','k','Tag','axesMagnitude');
        set(handles.axesReal,'Color','k','Tag','axesReal');
        set(handles.axesPhase,'Color','k','Tag','axesPhase');
        set(handles.axesImag,'Color','k','Tag','axesImag');
        
        if isempty(td)
            td = maxTimestampSec+ts;
        end
        
        % save the data
        setappdata(guiHandle,'timestampsSecs',td);
        setappdata(guiHandle,'digTimeDomainData',D);
        
    end
catch me
    % intentionally left blank
    fprintf('%s:%s\n',me.identifier,me.message);
end

% --- Executes during object deletion, before destroying properties.
function figRadix2FFT_DeleteFcn(hObject, ~, handles)
% hObject    handle to figRadix2FFT (see GCBO)
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


% --- Executes when user attempts to close figRadix2FFT.
function figRadix2FFT_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figRadix2FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
selection = questdlg('Close Radix2FFT?',...
                     'Close Request Function',...
                     'Yes','No','Yes');
switch selection,
   case 'Yes',
     delete(hObject);
   case 'No'
     return
end


% --------------------------------------------------------------------
function uitoggletoolNoise_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hfig      = handles.figRadix2FFT;
setappdata(hfig,'addNoise',~getappdata(hfig,'addNoise'));

% --- Executes on button press in checkboxAnalogData.
function checkboxAnalogData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAnalogData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAnalogData
guiHandle = handles.figRadix2FFT;
displayOn = getappdata(guiHandle,'displayOn');
displayOn(1) = ~displayOn(1);
setappdata(guiHandle,'displayOn',displayOn);


% --- Executes on button press in checkboxDigitalData.
function checkboxDigitalData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDigitalData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDigitalData
guiHandle = handles.figRadix2FFT;
displayOn = getappdata(guiHandle,'displayOn');
displayOn(2) = ~displayOn(2);
setappdata(guiHandle,'displayOn',displayOn);

% --- Executes on button press in checkboxFFTData.
function checkboxFFTData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFFTData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFFTData
guiHandle = handles.figRadix2FFT;
displayOn = getappdata(guiHandle,'displayOn');
displayOn(3) = ~displayOn(3);
setappdata(guiHandle,'displayOn',displayOn);

% --- Executes on button press in checkboxPhaseData.
function checkboxPhaseData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxPhaseData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxPhaseData
guiHandle = handles.figRadix2FFT;
displayOn = getappdata(guiHandle,'displayOn');
displayOn(5) = ~displayOn(5);
setappdata(guiHandle,'displayOn',displayOn);

% --- Executes on button press in checkboxRealData.
function checkboxRealData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRealData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxRealData
guiHandle = handles.figRadix2FFT;
displayOn = getappdata(guiHandle,'displayOn');
displayOn(4) = ~displayOn(4);
setappdata(guiHandle,'displayOn',displayOn);

% --- Executes on button press in checkboxImagData.
function checkboxImagData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxImagData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxImagData
guiHandle = handles.figRadix2FFT;
displayOn = getappdata(guiHandle,'displayOn');
displayOn(6) = ~displayOn(6);
setappdata(guiHandle,'displayOn',displayOn);

function resetAppData(hObject)
N = getappdata(hObject,'blockSize');
Fs = getappdata(hObject,'sampleRateHz');
setappdata(hObject,'timestampsSecs',[]);
setappdata(hObject,'digTimeDomainData',[]);
setappdata(hObject,'spectrogramImgQntzd',uint8(zeros(N,512,3)));
setappdata(hObject,'spectrogramImgDb',zeros(N,512));
setappdata(hObject,'spectrogramImgCol',1);

minDb = -90;
maxDb = 100;
dbLevels = linspace(minDb,maxDb,255);

setappdata(hObject,'spectrogramDbLevels',dbLevels);
setappdata(hObject,'spectrogramClrMap',uint8(255*colormap(hot(255))));

handles = guidata(hObject);

set(handles.sliderMinDb,'Min',-90,'Max',0,'Value',-90);
set(handles.sliderMaxDb,'Min',0,'Max',90,'Value',90);
set(handles.sliderSpecMaxFreqHz,'Min',8,'Max',Fs,'Value',N);


% --- Executes on slider movement.
function sliderMinDb_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMinDb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderMinDb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMinDb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderMaxDb_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMaxDb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderMaxDb_CreateFcn(hObject, ~, ~)
% hObject    handle to sliderMaxDb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderMaxDbContValCallback(guiHandle,~)

handles = guidata(guiHandle);
maxDb   = get(handles.sliderMaxDb,'Value');
minDb   = get(handles.sliderMinDb,'Value');

dbLevels = linspace(minDb,maxDb,255);
 
specImgQntzd    = getappdata(guiHandle,'spectrogramImgQntzd');
specImgDb       = getappdata(guiHandle,'spectrogramImgDb');
specImgCol      = getappdata(guiHandle,'spectrogramImgCol');
specImgClrMap   = getappdata(guiHandle,'spectrogramClrMap');

specImgDbLevels = dbLevels;

% end column
specImgCol = min(specImgCol,size(specImgDb,2));

N = size(specImgDb,1);

% quantize the data
for u=1:specImgCol

    for v=1:N
        idx = find(specImgDbLevels>specImgDb(v,u),1);
        if isempty(idx)
           specImgQntzd(v,u,:) = specImgClrMap(end,:);
        else
           specImgQntzd(v,u,:) = specImgClrMap(idx,:);
        end
    end
end

setappdata(guiHandle,'spectrogramImgQntzd',specImgQntzd);
setappdata(guiHandle,'spectrogramDbLevels',specImgDbLevels);

yTicks      = get(handles.axesSpectrogram,'YTick');
yTickLabels = get(handles.axesSpectrogram,'YTickLabel');
xTicks      = get(handles.axesSpectrogram,'XTick');
xTickLabels = get(handles.axesSpectrogram,'XTickLabel');

image(specImgQntzd,'Parent',handles.axesSpectrogram);
set(handles.axesSpectrogram,...
    'YTick',yTicks,...
    'YTickLabel',yTickLabels, ...        
    'XTick',xTicks,...
    'XTickLabel',xTickLabels, ...                      
    'Tag','axesSpectrogram');

function sliderSpecMaxFreqContValCallback(guiHandle,~)

handles      = guidata(guiHandle);
fs           = getappdata(guiHandle,'sampleRateHz');
specFreqMax  = floor(get(handles.sliderSpecMaxFreqHz,'Value'));
specImgQntzd = getappdata(guiHandle,'spectrogramImgQntzd');
xTicks       = get(handles.axesSpectrogram,'XTick');
xTickLabels  = get(handles.axesSpectrogram,'XTickLabel');
image(specImgQntzd(1:specFreqMax,:,:),'Parent',handles.axesSpectrogram);
set(handles.axesSpectrogram,...
    'YTick',linspace(1,specFreqMax,8),...
    'YTickLabel',num2str(floor(linspace(1,specFreqMax,8)')), ...          
    'XTick',xTicks,...
    'XTickLabel',xTickLabels, ...                     
    'Tag','axesSpectrogram');


% --- Executes during object creation, after setting all properties.
function sliderSpecMaxFreqHz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSpecMaxFreqHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
