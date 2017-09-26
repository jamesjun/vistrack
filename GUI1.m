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

% Last Modified by GUIDE v2.5 19-Mar-2014 15:24:51

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
csSettings = importdata('settings.m', '\n');
set(handles.editSettings, 'String', csSettings);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


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

[FileName,PathName,FilterIndex] = uigetfile('*.wmv;*.avi;*.mpg;*.mp4', ...
        'Select video file',get(handles.edit1, 'String'));
if FilterIndex
    try
        handles.vidFname = fullfile(PathName, FileName);
        set(handles.edit1, 'String', handles.vidFname);
        h = msgbox('Loading... (this will close automatically)');
        handles.vidobj = VideoReader(handles.vidFname);
        handles.vidobj
        try close(h); catch, end;   
        set(handles.btnSync, 'Enable', 'on');
        set(handles.btnBackground, 'Enable', 'off');
        set(handles.btnTrack, 'Enable', 'off');
        set(handles.btnPreview, 'Enable', 'off');
        set(handles.btnSave, 'Enable', 'off');
        set(handles.panelPlot, 'Visible', 'off');            
        msgstr = 'Video';
        % set the ADC file and ADC timestamp paths
        [~, fname, ~] = fileparts(FileName);
        handles.ADCfile = [PathName, fname, '_Rs.mat'];
        handles.ADCfileTs = [PathName, fname, '_Ts.mat'];
        set(handles.editADCfile, 'String', handles.ADCfile);
        set(handles.editADCfileTs, 'String', handles.ADCfileTs);        
        try handles.ADC = load(handles.ADCfile);
            msgstr = [msgstr, ', ADC_Rs'];
        catch, errordlg('ADC_Rs load error'); end
        try handles.ADCTS = load(handles.ADCfileTs);
            msgstr = [msgstr, ', ADC_Ts'];
        catch, errordlg('ADC_Ts load error'); end        
        
        guidata(hObject, handles);
        msgbox([msgstr ' file(s) loaded']);  
    catch
        errordlg(lasterr);
    end
end



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
FLIM0 = handles.FLIM0;

% Get time range from spike2
if ~exist('TLIM')
    try
        ADCTC = load(get(handles.editADCfileTs, 'String'));
        prefix = getSpike2Prefix(ADCTC);
        chTEXT = getfield(ADCTC, sprintf('%s_Ch%d', prefix, ADC_CH_TEXT));
        [EODR, TEOD, chName] = getSpike2Chan(handles.ADC, ADC_CH_EODR);
        AMPL = getSpike2Chan(handles.ADC, ADC_CH_AMPL);
        hfig = figure; AX = [];
        subplot 212; plot(TEOD, AMPL); AX(2) = gca; grid on;
        subplot 211; plot(TEOD, EODR); AX(1) = gca; grid on;        
        linkaxes(AX, 'x');
        hold on;                
        xlabel('Time (s)'); ylabel('EOD Rate (Hz)'); axis tight;                
        title({'Set time range and double-click', ...
               'r: GATE OPEN, m: ENTERED ARENA; g: GATE CLOSE', ...
               'c: FOUND FOOD; b: LIGHT BLINK; k: Default'}); 
        set(hfig, 'Name', handles.vidFname);
        TLIM = [nan, nan];
        for i=1:numel(chTEXT.times)            
            if ~isempty(strfind(chTEXT.text(i,:), 'GATE_OPEN'));                
                color = '-r';
            elseif ~isempty(strfind(chTEXT.text(i,:), 'ENTERED_ARENA'));
                color = '-m';                
            elseif ~isempty(strfind(chTEXT.text(i,:), 'GATE_CLOSE'));
                color = '-g';
            elseif ~isempty(strfind(chTEXT.text(i,:), 'FOUND_FOOD'));
                color = '-c';            
            elseif ~isempty(strfind(chTEXT.text(i,:), 'LIGHT_BLINK'));
                color = '-b';
            else
                color = '-k';
            end
            plot(chTEXT.times(i)*[1 1], get(gca, 'YLim'), color);
        end        
        gcax = get(gca, 'XLim'); 
        gcay = get(gca, 'YLim');
        h = imrect(gca, [gcax(1) gcay(1) diff(gcax) diff(gcay)]);
        hpos = wait(h);
        TLIM(1) = hpos(1);
        TLIM(2) = sum(hpos([1 3]));
        
        fprintf('TLIM: ');
        disp(TLIM(:)');
        try close(hfig), catch, end;
    catch
        errordlg('Specify TLIM = [First, Last]; in the Settings');
        disp(lasterr);
        handles.ADC
        return;
    end
end

% Set time range to track
FLIM1 = round(interp1(handles.TLIM0, FLIM0, TLIM, 'linear', 'extrap'));
FLIM1(1) = max(FLIM1(1), 1);
FLIMa = FLIM1(1) + [-149,150];   FLIMa(1) = max(FLIMa(1), 1);
FLIMb = FLIM1(2) + [-149,150];   FLIMb(1) = max(FLIMb(1), 1);

% Refine the first and last frames to track
try
    % first frame to track
    h=msgbox('Loading... (this will close automatically)');
    trImg = read(handles.vidobj, FLIMa);
    trImg = trImg(:,:,1,:);
    try close(h); catch, end;    
    implay(trImg);       
    uiwait(msgbox('Find the first frame to track and background, and close the movie'));
    answer = inputdlg({'First frame', 'Background frame'}, 'Get frames', 1, {'150', '300'});
    firstFrame = str2num(answer{1});
    img1 = trImg(:, :, 1, firstFrame);
    FLIM(1) = firstFrame + FLIM1(1) - 150;    
    bgFrame = str2num(answer{2});
    img00 = trImg(:, :, 1, bgFrame);    
    
    % last frame to track
    h=msgbox('Loading... (this will close automatically)');
    trImg = read(handles.vidobj, FLIMb);
    trImg = trImg(:,:,1,:);
    try close(h); catch, end;    
    implay(trImg);
    uiwait(msgbox('Find the last frame to track, and close the movie'));
    ans = inputdlg('Frame Number', 'Last frame', 1, {'150'});
    FLIM(2) = str2num(ans{1}) + FLIM1(2) - 150;

    % camera time unit
    TC = interp1(FLIM0, handles.TLIM0, FLIM(1):FLIM(2), 'linear');

    % Create background
    [img00, MASK, xy_init, vec0, xy0] = makeBackground(img1, img00);

catch
    disp(lasterr);
    return;
end

% Create a background
% ans = inputdlg({'Second image frame #'}, 'Background', 1, {sprintf('%0.0f', FLIM(2))});
% Frame2 = str2double(ans{1});
% h=msgbox('Loading... (this will close automatically)');
% img1 = read(handles.vidobj, Frame1); img1=img1(:,:,1);
% img2 = read(handles.vidobj, Frame2); img2=img2(:,:,1);
% try close(h); catch, end;
% [img00, MASK, xy_init, vec0, xy0] = makeBackground(img1, img2);

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
handles.fBWprev = TRACK_SHOW;
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

% Save
btnSave_Callback(hObject, eventdata, handles);


% --- Executes on button press in btnSync.
function btnSync_Callback(hObject, eventdata, handles)
% hObject    handle to btnSync (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;
vidobj = handles.vidobj;
if ~exist('FPS_cam')
    FPS_cam = get(vidobj, 'FrameRate');
end
if ~exist('TLIM0')
    try
        %load from spike2 TC channel
        ADCTC = load(get(handles.editADCfileTs, 'String'));
        prefix = getSpike2Prefix(ADCTC);
        chTCAM = getfield(ADCTC, sprintf('%s_Ch%d', prefix, ADC_CH_TCAM));
        chTCAM = chTCAM.times;
        TLIM0 = chTCAM([1 end]);
        fprintf('TLIM0: ');
        disp(TLIM0(:)');
    catch
        disp(lasterr);
        errordlg('Specify TLIM0 = [First, Last]; in the Settings');
        return;
    end
end

FLIM1 = [1, 300];
try
    h=msgbox('Loading... (this will close automatically)');    
    if ~isempty(vidobj.NumberOfFrames)
        FLIM1(2) = min(FLIM1(2), vidobj.NumberOfFrames);
    end    
    trImg = read(vidobj, FLIM1);
    trImg = trImg(:,:,1,:);    
    try close(h); catch, end;
    implay(trImg);
    uiwait(msgbox({'Find the first brightest blink, and close the movie', handles.vidFname}));
    ans = inputdlg('Frame Number', 'First frame',1,{'164'});
    FLIM0(1) = str2num(ans{1});

    FLIM0(2) = round(FLIM0(1) + FPS_cam * diff(TLIM0));
    FLIM1 = FLIM0(2) + [-100, 100]; 
    if ~isempty(vidobj.NumberOfFrames)
%         if FLIM1(1) < vidobj.NumberOfFrames
            FLIM1(2) = min(FLIM1(2), vidobj.NumberOfFrames);
%         end
    end 
    
    h=msgbox('Loading... (this will close automatically)');
    trImg = read(vidobj,FLIM1);
    trImg = trImg(:,:,1,:);    
    try close(h); catch, end;            
    implay(trImg);
    
    uiwait(msgbox({'Find the first brightest blink, and close the movie', handles.vidFname}));
    ans = inputdlg('Frame Number', 'Last frame',1,{'102'});
    temp = str2num(ans{1});
    FLIM0(2) = FLIM0(2)-100+temp-1;
catch
    disp(lasterr);
    disp(FLIM1);
    errordlg('Check the TLIM0 setting (ADC time range)');
    return;
end


% TLIM = interp1(FLIM0, TLIM0, FLIM([1 end]), 'linear');
FPS = diff(FLIM0)/diff(TLIM0);
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
img = handles.img1(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2));
img0c = img0(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2));
% img00c = handles.img00(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2));

dimg = uint8(img0c - img);

% absimg = imabsdiff(handles.img00, handles.img1);
% absimg(~handles.MASK) = 0;
% absimg = absimg(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2));
% figure; subplot 121; imagesc(absimg); title('absdiff');
% subplot 122; imagesc(dimg); title('diff');
% return;

BW = imdilate(bwmorph((dimg > IM_THRESH), 'clean', inf), SE);
BW = imfill(BW, 'holes');
BW = imclearborder(BW);

% dimg1 = dimg;
% dimg1(bwperim(BW)) = 256;
% imgabs = im(handles.img00, handles.img1);
% imgabs(~handles.MASK) = 0;
% figure; imagesc(imgabs);
% figure; imagesc(dimg1);
% return;

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
figure; 
subplot 221; 
imagesc(img, INTENSITY_LIM);  
axis equal; axis tight;
set(gca, {'XTick', 'YTick'}, {[],[]});
title('1. Original image');

subplot 222; 
imagesc(dimg); 
axis equal; axis tight; 
set(gca, {'XTick', 'YTick'}, {[], []});
title('2. Difference image');

subplot 223; imshow(BW);
title(sprintf('3. Binary image (Area=%d)', AreaTarget));    

subplot 224;  
BW1 = bwperim(BW);
img4 = img; img4(BW1)=255;     
% imshow(img4); 
imagesc(img4, INTENSITY_LIM);  
axis equal; axis tight;
set(gca, {'XTick', 'YTick'}, {[],[]});
title('4. Superimposed');
colormap gray;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

writeText('settings.m', get(handles.editSettings, 'String'));
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
writeText('settings.m', get(handles.editSettings, 'String'));
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

figure; imshow(gray2rgb(handles.img0)); title('Posture trajectory (red: more recent, circle: head)');
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


% --- Executes on button press in btnPlot4.
function btnPlot4_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlot4 (see GCBO)
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
figure; imagesc(handles.img0); 
set(gca, {'XTick', 'YTick'}, {[],[]}); axis equal; axis tight;
hold on; 
title(sprintf('EOD (%s) at the head trajectory (red: higher rate)', chName));

for i=1:numel(viColorRate)
    plotChevron([X2i(i), X3i(i)], [Y2i(i), Y3i(i)], mrColor(viColorRate(i),:), 90, .3);
%     plot(Xi(i), Yi(i), '.', 'color', mrColor(viColorRate(i),:));
end
colormap gray;


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
RLIM = [quantile(EODR1, .001), quantile(EODR1, .999)];
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


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;
h = msgbox('Calculating... (this will close automatically)');
%track head
[VISITCNT, TIMECNT] = calcVisitDensity(handles.img0, handles.TC, handles.XC(:,2), handles.YC(:,2), TRAJ_NFILT);
try close(h); catch, end;   

[~, exprID, ~] = fileparts(handles.vidFname);
figure; 
subplot 121; 
imshow(rgbmix(handles.img0, imgray2rgb((TIMECNT))));   
title(['Time density map: ', exprID]);

subplot 122; 
imshow(rgbmix(handles.img0, imgray2rgb((VISITCNT))));  
title(['Visit density map: ', exprID]);    

%update
handles.TIMECNT = TIMECNT;
handles.VISITCNT = VISITCNT;
guidata(hObject, handles);


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

[FileName,PathName,FilterIndex] = uigetfile('*_Track.mat','Select *_Track.mat file',get(handles.editResultFile, 'String'));
if FilterIndex
    try
        resultFile = fullfile(PathName, FileName);
        set(handles.editResultFile, 'String', resultFile);
        h = msgbox('Loading... (this will close automatically)');
        S = load(resultFile);
        try close(h); catch, end;        
        
        % Apply settings
        set(handles.editSettings, 'String', S.csSettings);
        csFields = {'TLIM0', 'FLIM0', 'FPS', ...
            'MASK' ,'xy_init' ,'vec0' ,'xy0' ,'TC' ,'TLIM' ,'FLIM' ,'img1' ,'img00', ...
            'SE' ,'thresh' ,'AreaTarget' ,'WINPOS' ,'img0', ...
            'XC' ,'YC' ,'AC' ,'xy_names' ,'ang_names' ,'csSettings', ...
            'ADC', 'ADCTS', ...
            'MOV', 'XC_off', 'YC_off', 'vidFname', 'ESAC'};
        for i=1:numel(csFields)
            eval(sprintf('handles.%s = S.%s;', csFields{i}, csFields{i}));
        end        
        
        set(handles.edit1, 'String', handles.vidFname);
        set(handles.editADCfile, 'String', [handles.vidFname(1:end-4), '_Rs.mat']);
        set(handles.editADCfileTs, 'String', [handles.vidFname(1:end-4), '_Ts.mat']);
        
        set(handles.btnSync, 'Enable', 'on');
        set(handles.btnBackground, 'Enable', 'on');
        set(handles.btnTrack, 'Enable', 'on');
        set(handles.btnPreview, 'Enable', 'on');
        set(handles.btnSave, 'Enable', 'on');
        set(handles.panelPlot, 'Visible', 'on'); 
        
        guidata(hObject, handles);
        msgbox('Tracking Result loaded');
    catch
        set(handles.editResultFile, 'String', '');
        errordlg(lasterr);
    end
end



% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ESAC = calcESAC(handles);

%save file
h = msgbox('Saving... (this will close automatically)');       
[pathstr, name, ext] = fileparts(handles.vidFname);
outfname = fullfile(pathstr, [name, '_Track.mat']);
eval(sprintf('save(''%s'', ''-struct'', ''handles'');', outfname));
try close(h); catch, end;
set(handles.editResultFile, 'String', outfname);
msgbox(sprintf('Output saved to %s', outfname));


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

FLIM = handles.FLIM;
MOV = handles.MOV;          
% figure; implay(MOV);return;

% figure; implay(MOV);
TC1 = handles.TC - handles.TC(1);
mrXC = bsxfun(@minus, handles.XC, handles.XC_off);
mrYC = bsxfun(@minus, handles.YC, handles.YC_off);
hfig = figure;
nF = size(MOV, 3);
tic1 = tic;
if ~exist('REPLAY_FLIM', 'var')
    REPLAY_FLIM = [1, nF];
end
iF=1;
try
    for iF=REPLAY_FLIM(1):REPLAY_STEP:REPLAY_FLIM(2)
        img1 = MOV(:,:,iF); 

        XC = mrXC(iF,:);
        YC = mrYC(iF,:);    

        % interpolated curve
        nxy = numel(XC);
        X1 = interp1(2:nxy, XC(2:end), 2:.1:nxy, 'spline');
        Y1 = interp1(2:nxy, YC(2:end), 2:.1:nxy, 'spline');

        clf(hfig);
        imshow(img1); hold on;
        figure(hfig);
        set(hfig, 'Name', handles.vidFname);
        plot(XC(2), YC(2), 'wo', XC(3:end), YC(3:end), 'ro',...
             X1, Y1, 'r-', XC(1), YC(1), 'g+'); %Mark the centroid
        title(sprintf('F1: %d; T1: %0.3f s', iF, TC1(iF)));
        drawnow;
        
        if exist('REPLAY_PAUSE', 'var')
            if REPLAY_PAUSE == 1
                pause;
            end
        end
    end
catch
    disp(lasterr);
    iF0 = iF + FLIM(1) - 1;
    tc0 = TC1(iF) + handles.TC(1);
    tc1 = TC1(iF);
    msgbox(sprintf('Closed at F0: %d, T0: %0.3f s; F1: %d, T1: %0.3f s\n', ...
            iF0, tc0, iF, tc1));
%     disp(lasterr);    
end
fprintf('Replay took %0.3f s; Realtime %0.3f s\n', toc(tic1), diff(TC1([1 end])));


% --- Executes on button press in btnFlip.
function btnFlip_Callback(hObject, eventdata, handles)
% hObject    handle to btnFlip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LOADSETTINGS;
ans = inputdlg('Frame #:', 'Flip orientation', 1, {''});
ans = str2double(ans);
if isnan(ans) || isempty(ans)
    return;
else
    viF = ans:numel(handles.TC);
end

%Flip X, Y
%xy_names = {'CoM', 'head', 'head-mid', 'mid', 'tail-mid', 'tail'};
handles.XC(viF, 2:6) = handles.XC(viF, 6:-1:2);
handles.YC(viF, 2:6) = handles.YC(viF, 6:-1:2);

%ang_names = {'CoM', 'head-mid', 'tail-mid', 'body-bend', 'tail-bend'};
AC(viF, :) = handles.AC(viF, :) + 180;
AC = mod(AC, 360);
AC(AC>180) = AC(AC>180) - 360;
handles.AC = AC;

%update
guidata(hObject, handles);
msgbox('Orientation Flipped');

% Save
btnSave_Callback(hObject, eventdata, handles);


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
