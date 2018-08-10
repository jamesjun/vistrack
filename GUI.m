function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 10-Aug-2018 11:32:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update settings window
% csSettings = importdata('settings_vistrack.m', '\n');
csSettings = file2cellstr_('settings_vistrack.m');
set(handles.editSettings, 'String', csSettings);

[vcVer, vcVer_date] = version_();
set(handles.textVer, 'String', sprintf('%s (%s) James Jun', vcVer, vcVer_date));
set(handles.btnUpdate, 'Enable', ifeq_(exist_dir_('.git'), 'on', 'off'));
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function csLines = file2cellstr_(vcFile)
% read text file to a cell string
fid = fopen(vcFile, 'r');
csLines = {};
while ~feof(fid), csLines{end+1} = fgetl(fid); end
fclose(fid);
csLines = csLines';


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnLoadVideo.
function btnLoadVideo_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles   in: {handles.edit1}
%           out: {vidfile, vidobj}

vidFname = get(handles.edit1, 'String');
if ~exist_file_(vidFname)
    [FileName,PathName,FilterIndex] = uigetfile('*.wmv;*.avi;*.mpg;*.mp4', ...
            'Select video file', vidFname);
    if ~FilterIndex, return; end
    vidFname = fullfile(PathName, FileName);
end
handles.vidFname = vidFname;
buttons_off_(handles);
clear_cache_();

% Set result field
vcFile_out = subsFileExt_(vidFname, '_Track.mat');
if exist_file_(vcFile_out)
    set(handles.editResultFile, 'String', vcFile_out);
    if ask_user_('Load previous tracking result?')
        try
            btnLoadPrev_Callback(handles.btnLoadPrev, eventdata, handles);
            set(handles.btnSync, 'Enable', 'off');
            set(handles.btnBackground, 'Enable', 'off');
            set(handles.btnTrack, 'Enable', 'off');
            set(handles.btnPreview, 'Enable', 'off');            
            set(handles.btnSave, 'Enable', 'on');
            set(handles.panelPlot, 'Visible', 'on');  
            return;
        catch
            msgbox('Loading tracking result failed');
        end
    end
else
    set(handles.editResultFile, 'String', '');
end

try      
    set(handles.edit1, 'String', handles.vidFname);
    h = msgbox('Loading... (this will close automatically)', 'modal');
%     handles.vidobj = VideoReader(handles.vidFname);
    handles.vidobj = vistrack('VideoReader', vidFname);
    fprintf('Video file info: %s\n', handles.vidFname);
    disp(handles.vidobj);
%     fprintf('Calculating number of frames...\n');
    fprintf('\t# video frames: %d\n', handles.vidobj.NumberOfFrames);
    close_(h);
    
    set(handles.btnSync, 'Enable', 'on');
    set(handles.btnBackground, 'Enable', 'off');
    set(handles.btnTrack, 'Enable', 'off');
    set(handles.btnPreview, 'Enable', 'off');
    set(handles.btnSave, 'Enable', 'off');
    set(handles.panelPlot, 'Visible', 'off');            
    msgstr = 'Video';
    % set the ADC file and ADC timestamp paths
%     [PathName, fname, ~] = fileparts(vidFname);
    
    % set timestamps if exists
    vcFile_Rs = subsFileExt_(vidFname, '_Rs.mat');
    if ~exist_file_(vcFile_Rs), vcFile_Rs = ''; end
    handles.ADCfile = vcFile_Rs;
    set(handles.editADCfile, 'String', vcFile_Rs);
    handles.ADC = try_load_(handles.ADCfile);

    vcFile_Ts = subsFileExt_(vidFname, '_Ts.mat');
    if ~exist_file_(vcFile_Ts), vcFile_Ts = ''; end
    handles.ADCfileTs = vcFile_Ts;
    set(handles.editADCfileTs, 'String', vcFile_Ts);   
    handles.ADCTS = try_load_(handles.ADCfileTs);
    
    if isempty(handles.ADC) || isempty(handles.ADCTS)
        set(handles.btnBackground, 'Enable', 'on');
        set(handles.btnSync, 'Enable', 'off');
    end
    
    guidata(hObject, handles);
    msgbox([msgstr ' file(s) loaded']);  
catch
    errordlg(lasterr);
end


function buttons_off_(handles)
set(handles.btnSync, 'Enable', 'off');
set(handles.btnBackground, 'Enable', 'off');
set(handles.btnTrack, 'Enable', 'off');
set(handles.btnPreview, 'Enable', 'off');
set(handles.btnSave, 'Enable', 'off');
set(handles.panelPlot, 'Visible', 'off');   


function editSettings_Callback(hObject, eventdata, handles)
% hObject    handle to editSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSettings as text
%        str2double(get(hObject,'String')) returns contents of editSettings as a double


% --- Executes during object creation, after setting all properties.
function editSettings_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    in: editSettings


% --- Executes on button press in btnBackground.
function btnBackground_Callback(hObject, eventdata, handles)
% hObject    handle to btnBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;

% Get time range from spike2  
vcFile_adc_ts = get(handles.editADCfileTs, 'String');
fSkipAdc = ~exist_file_(vcFile_adc_ts);
% if fSkipAdc
try
    [FLIM, TC, img1, img00] = mov_flim_(handles.vidobj);
    [img00, MASK, xy_init, vec0, xy0] = makeBackground(img1, img00);
    try
        TC = vistrack('cam2adc-sync', handles.S_sync, FLIM(1):FLIM(2));
%         TC = interp1(handles.FLIM0, handles.TLIM0, FLIM(1):FLIM(2), 'linear');
    catch
        ;
    end
catch
    return;
end

% Update handles structure
handles.MASK = MASK;
handles.xy_init = xy_init;
handles.vec0 = vec0;
handles.xy0 = xy0;
handles.TC = TC; 
handles.TLIM = TC([1, end]);
handles.FLIM = FLIM;
handles.img1 = img1;
handles.img00 = img00;
guidata(hObject, handles);
set(handles.btnPreview, 'Enable', 'on');


% --- Executes on button press in btnTrack.
function btnTrack_Callback(hObject, eventdata, handles)
% hObject    handle to btnTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;
handles.SE = strel('disk', BW_SMOOTH,0); %use disk size 3 for 640 480
handles.img0 = handles.img00 * (1-IM_MINCONTRAST); %make it darker
handles.fShow = TRACK_SHOW;
handles.ThreshLim = ThreshLim;
handles.fBWprev = 1;
TC = handles.TC;

try
    [XC, YC, AC, Area, S1, MOV, XC_off, YC_off] = trackFish(handles, handles.FLIM);
catch
    disp(lasterr)
    errordlg('Cancelled by user');
    return;
end

% Update figure handle
handles.XC = XC;
handles.YC = YC;
handles.AC = AC;
handles.MOV = MOV;
handles.XC_off = XC_off;
handles.YC_off = YC_off;
handles.xy_names = S1.xy_names;
handles.ang_names = S1.ang_names;
handles.csSettings = get(handles.editSettings, 'String');
set(handles.panelPlot, 'Visible', 'on');
set(handles.btnSave, 'Enable', 'on');
guidata(hObject, handles);
vistrack('trial-save', handles);


% --- Executes on button press in btnSync.
function btnSync_Callback(hObject, eventdata, handles)
% hObject    handle to btnSync (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = vistrack('trial-sync', handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in btnSync.
function btnSync_Callback1(hObject, eventdata, handles)

% S_sync = calc_sync_(handles);
LOADSETTINGS;
vidobj = handles.vidobj;
% mov1 = vistrack('loadvid-preview', handles.vidobj);

try
    handles.xyLED = xyLED;
catch
    handles.xyLED = round([vidobj.height, vidobj.width]/2);
end
if ~exist('FPS0', 'var'), FPS0 = get(vidobj, 'FrameRate'); end
ADCTC = load(get(handles.editADCfileTs, 'String'));
disp_adc_title_(ADCTC);
if ~exist('TLIM0', 'var')
    try
        %load from spike2 TC channel  
        prefix = getSpike2Prefix(ADCTC);
        chTCAM = getfield(ADCTC, sprintf('%s_Ch%d', prefix, ADC_CH_TCAM));
        chTCAM = chTCAM.times;
        if ~exist('SYNC_PERIOD', 'var'), SYNC_PERIOD = median(diff(chTCAM)); end
        TLIM0 = chTCAM([1 end]);
        fprintf('TLIM0: ');
        disp(TLIM0(:)');
    catch
        disp(lasterr);
        errordlg('Specify TLIM0 = [First, Last]; in the Settings');
        return;
    end
end
%     vistrack('checkLedPos', handles.vidobj, handles.xyLED);

try 
    if ~exist('FLIM0', 'var')
        FLIM0 = [];
        [FLIM0(1), handles.xyLed] = detectBlink(handles, 'first');
        FLIM0(2) = detectBlink(handles, 'last'); 
    end
    if ~exist('SYNC_FIRST', 'var'), SYNC_FIRST = 0; end
    FPS = diff(FLIM0) / diff(TLIM0);
    SYNC_SKIP0 = (diff(FLIM0)/FPS0 - diff(TLIM0))/SYNC_PERIOD;
    fprintf('TLIM: [%f, %f], SYNC_SKIP: %f\n', TLIM0(1), TLIM0(2), SYNC_SKIP0);
    SYNC_SKIP = round(SYNC_SKIP0);
    if SYNC_FIRST
        TLIM0(2) = TLIM0(2) + SYNC_SKIP*SYNC_PERIOD;
    else
        TLIM0(1) = TLIM0(1) - SYNC_SKIP*SYNC_PERIOD;
    end
    FPS = diff(FLIM0) / diff(TLIM0);
    fCheckSync = 0;
    
    %check if correct
    if fCheckSync
        [viBlink, vrBlinkInt] = findBlink(handles.vidobj, SYNC_PERIOD, FLIM0(1), FPS);    
        chTCAM = getfield(ADCTC, sprintf('%s_Ch%d', prefix, ADC_CH_TCAM));
        vtBlink = chTCAM.times;    
        n = min(numel(viBlink), numel(vtBlink));   
        viBlink2 = round(interp1(TLIM0, FLIM0, TLIM0(1):SYNC_PERIOD:TLIM0(end)));

        hFig = figure; 
        set(gcf, 'OuterPosition', get(0, 'ScreenSize')); drawnow;
        AX = [];
        subplot 311; plot(vtBlink(1:n), viBlink(1:n), 'b.', vtBlink(1:n), viBlink(1:n), 'bo', vtBlink(1:n), viBlink2, 'r.-'); 
        AX(1) = gca; title('Blue: Observed, Red: Predicted');
        xlabel('ADC Time (s)'); ylabel('Frame #');

        subplot 312; plot(vtBlink(1:n), viBlink(1:n) - viBlink2); ylabel('Frame# difference');
        AX(2) = gca;   

        subplot 313; plot(vtBlink(1:n), vrBlinkInt(1:n)); ylabel('Brightness');
        set(gca, 'YLim', [1 2]);
        AX(3) = gca;       

        linkaxes(AX, 'x');   
    end
catch
    disp(lasterr);
    errordlg('Check the TLIM0 setting (ADC time range)');
    return;
end
TLIM0 = round(TLIM0);
str = sprintf('FPS = %0.6f Hz, TLIM0=[%d, %d], FLIM0=[%d, %d]\n', FPS, TLIM0(1), TLIM0(2), FLIM0(1), FLIM0(2));
msgbox(str);
disp(str);

% Update handles structure
handles.TLIM0 = TLIM0;
handles.FLIM0 = FLIM0;
handles.FPS = FPS;
guidata(hObject, handles);
set(handles.btnBackground, 'Enable', 'on');


% --- Executes on button press in btnPreview.
function btnPreview_Callback(hObject, eventdata, handles)
% hObject    handle to btnPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;

% Obrain mask and initial location
img0 = handles.img00 * (1-IM_MINCONTRAST); %make it darker

% Initial
SE = strel('disk', BW_SMOOTH,0); %use disk size 3 for 640 480
[WINPOS, ~] = getBoundingBoxPos(handles.xy_init, size(img0), winlen*[1 1]);
img = handles.img1(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2), :);
img0c = img0(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2), :);
dimg = uint8_diff(img0c, img, 0);
BW = imdilate(bwmorph((dimg > IM_THRESH), 'clean', inf), SE);
BW = imfill(BW, 'holes');
BW = imclearborder(BW);

regions = regionprops(BW, {'Area', 'Centroid', 'Orientation', 'FilledImage', 'BoundingBox'});
if numel(regions)>1
    L = bwlabel(BW, 8);
    [iRegion] = region_nearest(regions, handles.xy0 - WINPOS([1,3]));
    regions = regions(iRegion);
    BW = L==iRegion;    
end
[BW, AreaTarget] = bwgetlargestblob(BW);

% Update
handles.SE = SE;
handles.thresh = IM_THRESH;
handles.AreaTarget = AreaTarget;
handles.WINPOS = WINPOS;
handles.img0 = img0;
guidata(hObject, handles);
set(handles.btnTrack, 'Enable', 'on');

% Preview images
hFig = figure;
set(hFig, 'OuterPosition', get(0, 'ScreenSize'), ...
    'Color', 'w', 'Name', handles.vidFname, 'NumberTitle', 'off');
drawnow;
subplot(2,2,1); 
imagesc(img, INTENSITY_LIM);  
axis equal; axis tight;
set(gca, {'XTick', 'YTick'}, {[],[]});
title('1. Original image');

subplot(2,2,2); 
imagesc(dimg); 
axis equal; axis tight; 
set(gca, {'XTick', 'YTick'}, {[], []});
title('2. Difference image');

subplot(2,2,3); 
imshow(BW);
title(sprintf('3. Binary image (Area=%d)', AreaTarget));    

subplot(2,2,4);  
BW1 = bwperim(BW);
img4 = img; img4(BW1)=255;     
% imshow(img4); 
imagesc(img4, INTENSITY_LIM);  
axis equal; axis tight;
set(gca, {'XTick', 'YTick'}, {[],[]});
title_('4. Superimposed (if incorrect, lower "IM_THRESH")');
colormap gray;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% do not auto save the settings unless manually saved
% writeText('settings.m', get(handles.editSettings, 'String'));
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
writeText('settings_vistrack.m', get(handles.editSettings, 'String'));
handles.csSettings = get(handles.editSettings, 'String');
guidata(hObject, handles);


function editADCfile_Callback(hObject, eventdata, handles)
% hObject    handle to editADCfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editADCfile as text
%        str2double(get(hObject,'String')) returns contents of editADCfile as a double


% --- Executes during object creation, after setting all properties.
function editADCfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editADCfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnLoadADC.
function btnLoadADC_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadADC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select ADC file',get(handles.editADCfile, 'String'));
if FilterIndex
    try
        handles.ADCfile = fullfile(PathName, FileName);
        set(handles.editADCfile, 'String', handles.ADCfile);
        h = msgbox('Loading... (this will close automatically)');
        handles.ADC = load(handles.ADCfile);
        try close(h); catch, end;        
        guidata(hObject, handles);
        msgbox('ADC File loaded');
    catch
        set(handles.editADCfile, 'String', '');
        errordlg(lasterr);
    end
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnPlot2.
function btnPlot2_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LOADSETTINGS;
XC = handles.XC;
YC = handles.YC;

hFig = figure; 
imshow((handles.img0)); title('Posture trajectory (red: more recent, circle: head)');
resize_figure(hFig, [0,0,.5,1]);
hold on;
nframes = size(XC,1);
nxy = size(XC,2);
mrColor = jet(nframes);
for iframe=1:TRAJ_STEP:nframes
    XI = interp1(2:nxy, XC(iframe,2:end), 2:.1:nxy, 'spline');
    YI = interp1(2:nxy, YC(iframe,2:end), 2:.1:nxy, 'spline');
    plot(XI, YI, 'color', mrColor(iframe,:));
    plot(XI(1), YI(1), 'o', 'color', mrColor(iframe,:));
end


% --- Executes on button press in btnEodRate.
function btnEodRate_Callback(hObject, eventdata, handles)
% hObject    handle to btnEodRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;
TC = handles.TC;
try
    ADC = handles.ADC;
    [EODR, TEOD, chName] = getSpike2Chan(ADC, ADC_CH_PLOT);
    fprintf('Loaded Spike2 %s (Ch%d)\n', chName, ADC_CH_PLOT);
catch
    disp(lasterr);
    errordlg('Load ADC file');
    return;
end


%smooth and interpolate position
TCi = TC(1):(1/EODR_SR):TC(end);
X2i = interp1(handles.TC, filtPos(handles.XC(:,2), TRAJ_NFILT), TCi, 'spline', 'extrap'); 
Y2i = interp1(handles.TC, filtPos(handles.YC(:,2), TRAJ_NFILT), TCi, 'spline', 'extrap'); 
X3i = interp1(handles.TC, filtPos(handles.XC(:,3), TRAJ_NFILT), TCi, 'spline', 'extrap'); 
Y3i = interp1(handles.TC, filtPos(handles.YC(:,3), TRAJ_NFILT), TCi, 'spline', 'extrap'); 

%convert rate to 0..255 color level at the camera time
R = interp1(TEOD, EODR, TCi);
R = (R-min(R))/(max(R)-min(R));
viColorRate = ceil(R * 256);
viColorRate(viColorRate<=0)=1;
mrColor = jet(256);
% figure; plot(handles.TC, viColorRate);

% Plot the EOD color representation
hFig = figure; 
imagesc(handles.img0); 
set(gca, {'XTick', 'YTick'}, {[],[]}); axis equal; axis tight;
hold on; 
title(sprintf('EOD (%s) at the head trajectory (red: higher rate)', chName));

for i=1:numel(viColorRate)
    plotChevron([X2i(i), X3i(i)], [Y2i(i), Y3i(i)], mrColor(viColorRate(i),:), 90, .3);
%     plot(Xi(i), Yi(i), '.', 'color', mrColor(viColorRate(i),:));
end
colormap gray;
resize_figure(hFig, [0 0 .5 1]);


% --- Executes on button press in btnPlot1.
function btnPlot1_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LOADSETTINGS;
XC = handles.XC;
YC = handles.YC;

% filter position
nf = TRAJ_NFILT;
XCf_h = filtPos(XC(:,2),nf);    YCf_h = filtPos(YC(:,2),nf);
XCf_hm = filtPos(XC(:,3),nf);   YCf_hm = filtPos(YC(:,3),nf);
XCf_m = filtPos(XC(:,4),nf);    YCf_m = filtPos(YC(:,4),nf);
XCf_tm = filtPos(XC(:,5),nf);    YCf_tm = filtPos(YC(:,5),nf);
XCf_t = filtPos(XC(:,6),nf);    YCf_t = filtPos(YC(:,6),nf);

figure; imshow(handles.img0); title('Trajectory of the rostral tip');
hold on; plot(XCf_h, YCf_h);
pause(.4);
hold on; comet(XCf_h, YCf_h, .1);


% --- Executes on button press in btnPlot3.
function btnPlot3_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlot3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnMovOut.
function btnMovOut_Callback(hObject, eventdata, handles)
% hObject    handle to btnMovOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load data
LOADSETTINGS;
TC = handles.TC;
XC = handles.XC;
YC = handles.YC;
try
    ADC = handles.ADC;
    [EODR, TEOD, chName] = getSpike2Chan(ADC, ADC_CH_PLOT);
    fprintf('Loaded Spike2 %s (Ch%d)\n', chName, ADC_CH_PLOT);
    
    % Resample EOD Rate to camera time    
    RC = interp1(TEOD, EODR, TC);
%     RC = filtfilt(ones(1,TRAJ_STEP), TRAJ_STEP, RC);
catch
    disp(lasterr)
    errordlg('Load ADC file');
    return;
end
% figure; plot(TEOD, EODR, 'r.', TC, RC, 'b-');

% Plot the EOD color representation
writerObj = VideoWriter(MOV_FILEOUT, 'MPEG-4');
set(writerObj, 'FrameRate', handles.FPS); %30x realtime
% set(writerObj, 'Quality', 90); %30x realtime
open(writerObj);

figure; title('Reposition'); pause;
subplot(4,1,1:3);
hfig = imshow(gray2rgb(handles.img0, INTENSITY_LIM));
% hfig = imagesc(handles.img0, INTENSITY_LIM); 
set(gca, {'XTick', 'YTick'}, {[],[]}); 
axis equal; axis tight; hold on; 
title('EOD rate at the head (red: higher rate)');

%plot locations
nframes = size(XC,1);
% mrColor = jet(nframes);
[mrColor, vrRateSrt, vrQuantSrt] = quantile2color(RC);

%colorbar
plotColorbar(size(handles.img0), vrRateSrt, vrQuantSrt);
EODR1 = EODR(TEOD > TLIM(1) & TEOD < TLIM(2));
RLIM = [quantile_(EODR1, .001), quantile_(EODR1, .999)];
htext = [];
vhChevron = [];
for iframe=1:nframes
    
    %------------------
    subplot(4,1,1:3);
    frameNum = iframe + handles.FLIM(1) - 1;
    mrImg = readFrame(handles.vidobj, frameNum);    
    mrImg(~handles.MASK) = 0;   
    mrImg = gray2rgb(mrImg, INTENSITY_LIM);
    set(hfig, 'cdata', mrImg);    

    try delete(htext); catch, end;
    htext(1) = text(10, 30, sprintf('EOD (%s): %0.1f Hz', chName, RC(iframe)), ...
        'FontSize', 12, 'Color', [1 1 1]);
    htext(2) = text(10, 75, sprintf('Time: %0.1f s', TC(iframe)), ...
        'FontSize', 12, 'Color', [1 1 1]);    
    htext(3) = text(10, 120, sprintf('Frame: %d', frameNum), ...
        'FontSize', 12, 'Color', [1 1 1]);    
    
    if mod(iframe, MOV_PLOTSTEP) == 0
        vhChevron(end+1) = plotChevron(XC(iframe, 2:3), YC(iframe, 2:3), mrColor(iframe,:), 90, .3);
        if numel(vhChevron) > MOV_PLOTLEN
            delete(vhChevron(1:end-MOV_PLOTLEN));
            vhChevron(1:end-MOV_PLOTLEN) = [];
        end
    end
    
    %------------------
    subplot(4,1,4);    
    hold off;
    plot(TEOD - TC(iframe), EODR, 'k.'); hold on;
    axis([MOV_TimeWin(1), MOV_TimeWin(2), RLIM(1), RLIM(2)]);
    plot([0 0], get(gca, 'YLim'), 'r-');
    grid on;
    xlabel('Time (sec)'); 
    ylabel(sprintf('EOD (%s) Hz', chName));  
    
    colormap jet;
%     drawnow;
    try
        writeVideo(writerObj, getframe(gcf));
    catch
        disp('Movie output cancelled by user');
        close(writerObj);
        return;
    end
end

close(writerObj);
msgbox(sprintf('File written to %s', MOV_FILEOUT));

% --- Executes on button press in btnSoundOut
function btnSoundOut_Callback(hObject, eventdata, handles)
% hObject    handle to btnSoundOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load data
LOADSETTINGS;
TC = handles.TC;
TLIM = [TC(1), TC(end)];

try
    S = handles.ADCTS;
    csFieldnames = fieldnames(S);
    S = getfield(S, csFieldnames{1});    
    Teod = S.times;
    Teod1 = Teod(Teod > TLIM(1) & Teod < TLIM(2));
    viEod1 = round((Teod1 - TLIM(1)) * WAV_Fs);
    tdur = diff(TLIM);
    ns = round(tdur * WAV_Fs);
    
    %make a binary vector
    mlBinvec = zeros(ns,1);
    mlBinvec(viEod1) = 1;
    wavwrite(mlBinvec, WAV_Fs, WAV_FILEOUT);
    
    msgbox(sprintf('File written to %s', WAV_FILEOUT));
catch
    disp(lasterr)
    errordlg('Load ADC Timestamp');
    return;
end


function editADCfileTs_Callback(hObject, eventdata, handles)
% hObject    handle to editADCfileTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editADCfileTs as text
%        str2double(get(hObject,'String')) returns contents of editADCfileTs as a double


% --- Executes during object creation, after setting all properties.
function editADCfileTs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editADCfileTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnLoadADCTS.
function btnLoadADCTS_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadADCTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select ADC Timestamp',get(handles.editADCfileTs, 'String'));
if FilterIndex
    try
        handles.ADCfileTs = fullfile(PathName, FileName);
        set(handles.editADCfileTs, 'String', handles.ADCfileTs);
        h = msgbox('Loading... (this will close automatically)');
        handles.ADCTS = load(handles.ADCfileTs);
        try close(h); catch, end;        
        guidata(hObject, handles);
        msgbox('ADC File loaded');
    catch
        set(handles.editADCfileTs, 'String', '');
        errordlg(lasterr);
    end
end


% --- Executes during object creation, after setting all properties.
function btnPreview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btnPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btnPlotDensity.
function btnPlotDensity_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlotDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% vistrack('trial-visitcount', handles);
vistrack('trial-visitcount', handles);

% LOADSETTINGS;
% h = msgbox('Calculating... (this will close automatically)');
% 
% % [mnVisit1, mnVisit] = calcVisitCount(vsTrialPool, img0, mlMask, nGrid);
% 
% %track head
% [VISITCNT, TIMECNT] = calcVisitDensity(handles.img0, handles.TC, handles.XC(:,2), handles.YC(:,2), TRAJ_NFILT);
% try close(h); catch, end;   
% 
% [~, exprID, ~] = fileparts(handles.vidFname);
% img0_adj = imadjust(handles.img0);
% figure; 
% subplot 121; 
% imshow(rgbmix(img0_adj, imgray2rgb(TIMECNT), [], 'transparent'));   
% title(['Time density map: ', exprID]);
% 
% subplot 122; 
% imshow(rgbmix(img0_adj, imgray2rgb(VISITCNT), [], 'transparent'));   
% title(['Visit density map: ', exprID]);    
% 
% %update
% handles.TIMECNT = TIMECNT;
% handles.VISITCNT = VISITCNT;
% guidata(hObject, handles);


function editResultFile_Callback(hObject, eventdata, handles)
% hObject    handle to editResultFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editResultFile as text
%        str2double(get(hObject,'String')) returns contents of editResultFile as a double


% --- Executes during object creation, after setting all properties.
function editResultFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editResultFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnLoadPrev.
function btnLoadPrev_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resultFile = get(handles.editResultFile, 'String');
if ~exist_file_(resultFile)
    [FileName,PathName,FilterIndex] = uigetfile('*_Track.mat','Select *_Track.mat file', resultFile);
    if ~FilterIndex, return; end
    resultFile = fullfile(PathName, FileName);
end
try
    %   resultFile = fullfile(PathName, FileName);
    set(handles.editResultFile, 'String', resultFile);
    h = msgbox('Loading... (this will close automatically)');
    warning off;
    S = load(resultFile);
    try close(h); catch, end;        

    % Apply settings
    set(handles.editSettings, 'String', S.csSettings);
    S_cfg = vistrack('load-cfg');
    handles = struct_merge_(handles, S, S_cfg.csFields);
    handles.vcFile_Track = resultFile;

    set(handles.edit1, 'String', handles.vidFname);
    set(handles.editADCfile, 'String', [handles.vidFname(1:end-4), '_Rs.mat']);
    set(handles.editADCfileTs, 'String', [handles.vidFname(1:end-4), '_Ts.mat']);

    set(handles.btnSync, 'Enable', 'off');
    set(handles.btnBackground, 'Enable', 'off');
    set(handles.btnTrack, 'Enable', 'off');
    set(handles.btnPreview, 'Enable', 'off');
    set(handles.btnSave, 'Enable', 'on');
    set(handles.panelPlot, 'Visible', 'on'); 

    guidata(hObject, handles);
    clear_cache_();
    msgbox('Tracking Result loaded');
catch
    set(handles.editResultFile, 'String', '');
    errordlg(lasterr);
end



% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trial-save', handles);

% --- Executes during object creation, after setting all properties.
function btnSave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btnReplay.
function btnReplay_Callback(hObject, eventdata, handles)
% hObject    handle to btnReplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning off;
LOADSETTINGS;

MOV = handles.MOV;          
TC1 = handles.TC - handles.TC(1);
mrXC = bsxfun(@minus, handles.XC, handles.XC_off);
mrYC = bsxfun(@minus, handles.YC, handles.YC_off);

%setup fig
timer1 = timer('Period', 1/15, 'ExecutionMode', 'fixedRate', 'TasksToExecute', inf);
hFig = figure('NumberTitle', 'off', 'Name', ...
    '[H]elp; (Sft)+[L/R/U/D]; SPACEBAR:Pause; [F]lip; [C]ut upto here; [G]oto Frame');
hImg = imshow(imadjust(MOV(:,:,1))); hold on;
iFrame = 1;
XC = mrXC(iFrame,:);
YC = mrYC(iFrame,:);
nxy = numel(XC);
X1 = interp1(2:nxy, XC(2:end), 2:.1:nxy, 'spline');
Y1 = interp1(2:nxy, YC(2:end), 2:.1:nxy, 'spline');
hPlot = plot(XC(2), YC(2), 'go', XC(3:end), YC(3:end), 'mo',...
     X1, Y1, 'm-', XC(1), YC(1), 'g+', 'LineWidth', 2); %Mark the centroid
hTitle = title(sprintf('F1: %d; T1: %0.3f s, Step: %d', ...
    iFrame, TC1(iFrame), REPLAY_STEP));
resize_figure(hFig, [0,0,.5,1]);

%setup timer
Stimer = struct('hImg', hImg, 'TC1', TC1, 'hFig', hFig, 'hPlot', hPlot, 'iFrame', iFrame, ...
    'hTitle', hTitle, 'REPLAY_STEP', REPLAY_STEP, 'hObject', hObject);
set(timer1, 'UserData', Stimer);
set(hFig, 'UserData', timer1, 'KeyPressFcn', @keyFcnPreview, ...
    'CloseRequestFcn', @closeFcnPreview);
set(timer1, 'TimerFcn', @timerFcnPreview);
start(timer1);


% --- Executes on button press in btnFlip.
function btnFlip_Callback(hObject, eventdata, handles)
% hObject    handle to btnFlip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LOADSETTINGS;
iFrame = str2double(inputdlg('Frame #:', 'Flip orientation', 1, {''}));
if isempty(iFrame), return ;end
GUI_FLIP;
% Save
% btnSave_Callback(hObject, eventdata, handles);


% --- Executes on button press in btnCustom.
function btnCustom_Callback(hObject, eventdata, handles)
% hObject    handle to btnCustom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_ESACPLOT;


% --- Executes on button press in btnESAC.
function btnESAC_Callback(hObject, eventdata, handles)
% hObject    handle to btnESAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_ESACCADES;


% --- Executes on button press in btnEODMovie.
function btnEODMovie_Callback(hObject, eventdata, handles)
% hObject    handle to btnEODMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_EODMOVIE;


function msgbox_(vcMsg)
fprintf('%s\n', vcMsg);
msgbox(vcMsg);


function disp_adc_title_(ADCTC)
S_adc = ADCTC;
fprintf('Spike2 ADC output:\n');
cellfun(@(vc)fprintf('\t%s: %s (%0.3f-%0.3fs)\n', ...
    vc, S_adc.(vc).title, S_adc.(vc).times(1), S_adc.(vc).times(end)), ...
        fieldnames(S_adc));
fprintf('\n');


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile)
if nargin<2, fVerbose = 0; end
if isempty(vcFile)
    flag = 0; 
else
    S_dir = dir(vcFile);
    if numel(S_dir) == 1
        flag = ~S_dir.isdir;
    else
        flag = 0;
    end
end
if fVerbose && ~flag
    fprintf(2, 'File does not exist: %s\n', vcFile);
end


%--------------------------------------------------------------------------
% 8/2/17 JJJ: added '.' if dir is empty
% 7/31/17 JJJ: Substitute file extension
function varargout = subsFileExt_(vcFile, varargin)
% Substitute the extension part of the file
% [out1, out2, ..] = subsFileExt_(filename, ext1, ext2, ...)

[vcDir_, vcFile_, ~] = fileparts(vcFile);
if isempty(vcDir_), vcDir_ = '.'; end
for i=1:numel(varargin)
    vcExt_ = varargin{i};    
    varargout{i} = [vcDir_, filesep(), vcFile_, vcExt_];
end


%--------------------------------------------------------------------------
function flag = ask_user_(vcMsg, fYes)
if nargin<2, fYes = 1; end
if fYes
    vcAns = questdlg(vcMsg, '', 'Yes', 'No', 'Yes');
else
    vcAns = questdlg(vcMsg, '', 'Yes', 'No', 'Yes');
end
flag = strcmp(vcAns, 'Yes');


%--------------------------------------------------------------------------
function S = try_load_(vcFile)
S = [];
if isempty(vcFile), return; end
try
    S = load(vcFile);
    fprintf('%s loaded\n', vcFile);
catch
    errordlg(sprintf('%s load error', vcFile)); 
end


% --- Executes on button press in btnUpdate.
function btnUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to btnUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('update');
msgbox('Update successful, Restart the app');


%--------------------------------------------------------------------------
% 9/29/17 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate] = version_()
if nargin<1, vcFile_prm = ''; end
% vcVer = 'v0.1.2';
% vcDate = '7/11/2018';
[vcVer, vcDate] = vistrack('version');
if nargout==0
    fprintf('%s (%s) installed\n', vcVer, vcDate);
end


% --- Executes on button press in btnLoadTrialSet.
function btnLoadTrialSet_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadTrialSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vcFile_trialset = get(handles.editTrialSet, 'String');
if ~exist_file_(vcFile_trialset)
    [FileName,PathName,FilterIndex] = uigetfile('*.trialset', ...
            'Select trialset', vcFile_trialset);
    if ~FilterIndex, return; end
    vcFile_trialset = fullfile(PathName, FileName);
else
    vcFile_trialset = fullpath_(vcFile_trialset);
end
set(handles.editTrialSet, 'String', vcFile_trialset);
set(handles.panelTrialSet, 'Visible', 'on');
edit(vcFile_trialset);
msgbox(sprintf('Loaded %s', vcFile_trialset), 'Loading Trialset');


% --- Executes on button press in btnEditTrialset.
function btnEditTrialset_Callback(hObject, eventdata, handles)
% hObject    handle to btnEditTrialset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit(get(handles.editTrialSet, 'String'));


% --- Executes on button press in btnLearningCurve.
function btnLearningCurve_Callback(hObject, eventdata, handles)
% hObject    handle to btnLearningCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-learningcurve', get_str_(handles.editTrialSet)); 
% [cvrPathLen, cvrDuration] = vistrack('measure_trials', get_str_(handles.editTrialSet)); 


% --- Executes on button press in btnSamplingDensity.
function btnSamplingDensity_Callback(hObject, eventdata, handles)
% hObject    handle to btnSamplingDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnESCAN.
function btnESCAN_Callback(hObject, eventdata, handles)
% hObject    handle to btnESCAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnBSCAN.
function btnBSCAN_Callback(hObject, eventdata, handles)
% hObject    handle to btnBSCAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnVisitDensity.
function btnVisitDensity_Callback(hObject, eventdata, handles)
% hObject    handle to btnVisitDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnProbeTrials.
function btnProbeTrials_Callback(hObject, eventdata, handles)
% hObject    handle to btnProbeTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-probe', get_str_(handles.editTrialSet));


% --- Executes on button press in btnListFiles.
function btnListFiles_Callback(hObject, eventdata, handles)
% hObject    handle to btnListFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-list', get_str_(handles.editTrialSet));


function vc = get_str_(hObj)
try
    vc = get(hObj, 'String');
catch
    vc = '';
end


% --- Executes on button press in btnBarPlots.
function btnBarPlots_Callback(hObject, eventdata, handles)
% hObject    handle to btnBarPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-barplots', get_str_(handles.editTrialSet));


% --- Executes on button press in btnExportCsv.
function btnExportCsv_Callback(hObject, eventdata, handles)
% hObject    handle to btnExportCsv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-exportcsv', get_str_(handles.editTrialSet));


% --- Executes on button press in btnTrialsetFps.
function btnTrialsetFps_Callback(hObject, eventdata, handles)
% hObject    handle to btnTrialsetFps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-checkfps', get_str_(handles.editTrialSet));


% --- Executes on button press in btnTrialsetCoordinates.
function btnTrialsetCoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to btnTrialsetCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-coordinates', get_str_(handles.editTrialSet));


% --- Executes on button press in btnFixFps_trialset.
function btnFixFps_trialset_Callback(hObject, eventdata, handles)
% hObject    handle to btnFixFps_trialset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('trialset-fixfps', get_str_(handles.editTrialSet));


% --- Executes on button press in pushbutton70.
function pushbutton70_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function [FLIM, TC, img1, img00] = mov_flim_(vidobj, nFrames_load)
% mov_flim_(): clear cache

persistent csAns1 csAns2 csAns3
if nargin==0, [csAns1, csAns2, csAns3] = deal([]); return; end %clear cache

if nargin<2, nFrames_load = 300; end % skip every 10 frames

% if isempty(nFrames_skip), nFrames_skip = 75; end
warning off;
nFrames = vidobj.NumberOfFrames;
% nFrames_skip = floor(nFrames / nFrames_load);
TC = linspace(0, vidobj.Duration, nFrames); % use sync function 
% viF = 1:nFrames_skip:nFrames;
viF = unique(round(linspace(1, nFrames, nFrames_load)));
% viF = 1:nFrames_load; % much faster to load
% if isempty(mov1)
mov1 = vistrack('loadvid-preview', vidobj);
%     mov1 = mov_shrink_(read_(vidobj, viF), 2); 
% end

% rough scan
implay(mov1(:,:,viF));
uiwait(msgbox('Find the first and last frame to track, and close the movie'));
if isempty(csAns1), csAns1 = {'1', sprintf('%d', numel(viF))}; end
csAns = inputdlg({'First frame', 'Last frame'}, 'Get frames', 1, csAns1);
csAns1 = csAns;
frame_first = viF(str2num(csAns{1}));
frame_last = viF(str2num(csAns{2}));
fprintf('#1: First frame: %s; Last frame: %s\n', csAns{1}, csAns{2});

% Find first frame to track
viF_first = trim_(frame_first + (-150:149), 1, nFrames);
% tmr = read_(vidobj, viF_first);
implay(mov1(:,:,viF_first));
% implay(mov_shrink_(tmr, 2));
uiwait(msgbox('Find the first frame to track and background, and close the movie'));
if isempty(csAns2), csAns2 = {'1', num2str(numel(viF_first))}; end
csAns = inputdlg({'First frame', 'Background frame'}, 'Get frames', 1, csAns2);
csAns2 = csAns;
frame_first = viF_first(str2num(csAns{1}));
frame_background = viF_first(str2num(csAns{2}));
img1 = read_(vidobj, frame_first);
img00 = read_(vidobj, frame_background);
fprintf('#2: First frame: %s; Background frame: %s\n', csAns{1}, csAns{2});

% Find last frame to track
viF_last = trim_(frame_last + (-150:149), 1, nFrames);
implay(mov1(:,:,viF_last));
uiwait(msgbox('Find the last frame to track, and close the movie'));
iFrame3 = ceil(numel(viF_last)/2);
if isempty(csAns3), csAns3 = {num2str(iFrame3)}; end
csAns = inputdlg({'Last frame'}, 'Get frames', 1, csAns3);
csAns3 = csAns;
frame_last = viF_last(str2num(csAns{1}));
fprintf('#3: Last frame: %s\n', csAns{1});

FLIM = [frame_first, frame_last];
TC = TC(frame_first:frame_last);


function close_(h)
try close(h); catch; end


function tmr = read_(vidobj, viF)
h=msgbox('Loading video... (this will close automatically)'); drawnow;
tmr = vid_read(vidobj, viF);
close_(h);


function vi=trim_(vi, a, b)
vi = vi(vi>=a & vi<=b);


% --- Executes on button press in btnImport_Track.
function btnImport_Track_Callback(hObject, eventdata, handles)
% hObject    handle to btnImport_Track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vcFile_trialset = get(handles.editTrialSet, 'String');
vistrack('trialset-import-track', vcFile_trialset);


% --- Executes on button press in btnOpenSheet.
function btnOpenSheet_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenSheet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vcFile_trialset = get(handles.editTrialSet, 'String');
vistrack('trialset-googlesheet', vcFile_trialset);


% --- Executes on button press in btnSummary.
function btnSummary_Callback(hObject, eventdata, handles)
% hObject    handle to btnSummary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
csMsg = vistrack('summary', handles);
disp_cs_(csMsg);
msgbox(csMsg);


% --- Executes on button press in btnExport.
function btnExport_Callback(hObject, eventdata, handles)
% hObject    handle to btnExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vistrack('export', handles);


%--------------------------------------------------------------------------
% display cell string
function disp_cs_(cs)
cellfun(@(s)fprintf('%s\n',s), cs);


%--------------------------------------------------------------------------
% Retreive full path of a file
function vcFile_full = fullpath_(vcFile)
[vcDir_, vcFile_, vcExt_] = fileparts(vcFile);
if isempty(vcDir_) 
    vcDir_ = pwd();
    vcFile_full = fullfile(vcDir_, vcFile);
else
    vcFile_full = vcFile;
end
if nargout>=2, S_dir = dir(vcFile_full); end


%--------------------------------------------------------------------------
% 7/20/2018 Copied from jrc3.m
% If field doesn't exist copy empty
function P = struct_merge_(P, P1, csNames)
% Merge second struct to first one
% P = struct_merge_(P, P_append)
% P = struct_merge_(P, P_append, var_list) : only update list of variable names
if isempty(P), P=P1; return; end % P can be empty
if isempty(P1), return; end
if nargin<3, csNames = fieldnames(P1); end
if ischar(csNames), csNames = {csNames}; end

for iField = 1:numel(csNames)
    vcName_ = csNames{iField};
    if isfield(P1, vcName_)
        P.(vcName_) = P1.(vcName_); 
    else
        P.(vcName_) = [];
    end
end


%--------------------------------------------------------------------------
% 7/24/2018 JJJ: Copied from jrc3.m 
function flag = exist_dir_(vcDir)
if isempty(vcDir)
    flag = 0;
else
    flag = exist(vcDir, 'dir') == 7;
end


%--------------------------------------------------------------------------
% 7/24/2018: Copied from jrc3.m
function flag = key_modifier_(event, vcKey)
% Check for shift, alt, ctrl press
try
    flag = any(strcmpi(event.Modifier, vcKey));
catch
    flag = 0;
end


%--------------------------------------------------------------------------
% 7/24/2018: Copied from jrc3.m
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end


%--------------------------------------------------------------------------
% 7/24/2018 JJJ: Clear persistent memories
function clear_cache_()
mov_flim_(); % clear cache;
vistrack('clear-cache');


%--------------------------------------------------------------------------
% 7/24/2018 JJJ: resize movie by a scale factor
function mov = mov_shrink_(mov, nSkip)
if nSkip==1, return; end
mov = mov(1:nSkip:end, 1:nSkip:end, :);


%--------------------------------------------------------------------------
function hTItle = title_(vc)
hTItle = title(vc, 'Interpreter', 'none');


% --- Executes on button press in btnFixSync.
function btnFixSync_Callback(hObject, eventdata, handles)
% hObject    handle to btnFixSync (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = vistrack('trial-fixsync', handles);
guidata(hObject, handles);


% --- Executes on button press in pushbutton78.
function pushbutton78_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton81.
function pushbutton81_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton82.
function pushbutton82_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
