function keyFcnPreview(hFig, event)

% ensure the figure is still valid
if ~ishandle(hFig), return; end
    
timer1 = get(hFig, 'UserData');
S = get(timer1, 'UserData');
fRunning =  strcmpi(timer1.Running, 'on');
handles = guidata(S.hObject);

nFrames = size(handles.MOV,3);
[~, vcDataID, ~] = fileparts(handles.vidFname);
[iFrame, hObject] = deal(S.iFrame, S.hObject);
switch lower(event.Key)
    case 'h' %help
        csHelp = { ...
            '----Playback----', 
            '[SPACE]: start and stop video', 
            '(Shift) + LEFT/RIGHT: backward/forward (Shift: 4x)', 
            '[UP/DOWN]: Speed up/down', 
            '[HOME/END]: Go to start/end',
            '[G]oto: Go to a specific frame}', 
            '----EDIT----',
            '[F]lip head/tail', 
            '[C]ut: Trim video up to the current frame'};
        msgbox(csHelp);
        
    case 'space'  % toggle start and stop  
        if fRunning, stop(timer1);
        else start(timer1); end
    
    case {'leftarrow', 'rightarrow', 'f', 'home', 'end', 'g'}
    %'f' for flip, 'left' for back, 'right' for forward
        if fRunning, stop(timer1); end
        switch event.Key
            case {'leftarrow', 'rightarrow'}
                nStep = S.REPLAY_STEP * ifeq_(key_modifier_(event, 'shift'), 4, 1);
                nStep = nStep * ifeq_(strcmpi(event.Key, 'leftarrow'), -1, 1);
                S.iFrame = min(max(1, S.iFrame + nStep), nFrames);
            case 'home', S.iFrame = 1;
            case 'end', S.iFrame = nFrames;
            case 'f', GUI_FLIP; %flip orientation
            case 'g' % go to frame
                S.iFrame = uigetnum(sprintf('Go to Frame (choose from 1-%d)', nFrames), S.iFrame);
                if isempty(S.iFrame), return; end
                if (S.iFrame < 1 || S.iFrame > nFrames), S.iFrame = nan; end
                if isnan(S.iFrame), msgbox('Cancelled'); return; end
        end  
        set(S.hImg, 'CData', imadjust(handles.MOV(:,:,S.iFrame)));
        
        %annotate iamge
        handles = guidata(S.hObject);
        mrXC = bsxfun(@minus, handles.XC, handles.XC_off);
        mrYC = bsxfun(@minus, handles.YC, handles.YC_off);
        XC = mrXC(S.iFrame,:);
        YC = mrYC(S.iFrame,:);
        nxy = numel(XC);
        X1 = interp1(2:nxy, XC(2:end), 2:.1:nxy, 'spline');
        Y1 = interp1(2:nxy, YC(2:end), 2:.1:nxy, 'spline');
        set(S.hPlot(1), 'XData', XC(2), 'YData', YC(2));
        set(S.hPlot(2), 'XData', XC(3:end), 'YData', YC(3:end));
        set(S.hPlot(3), 'XData', X1, 'YData', Y1);
        set(S.hPlot(4), 'XData', XC(1), 'YData', YC(1));
                
    case 'uparrow' %speed up
        S.REPLAY_STEP = min(S.REPLAY_STEP+1, 30);
        
    
    case 'downarrow' %speed up
        S.REPLAY_STEP = max(S.REPLAY_STEP-1, 1);  
        
    case 'c' % cut up to here
        if fRunning, stop(timer1); end
        GUI_CUT;
        
end %switch

set(S.hTitle, 'String', ...
    sprintf('F1: %d; T1: %0.3f s, Step: %d (%s)', ...
    S.iFrame, S.TC1(S.iFrame), S.REPLAY_STEP, vcDataID));

set(timer1, 'UserData', S);
end %func


%--------------------------------------------------------------------------
% 7/24/2018: Copied from jrc3.m
function flag = key_modifier_(event, vcKey)
% Check for shift, alt, ctrl press
try
    flag = any(strcmpi(event.Modifier, vcKey));
catch
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
% 7/24/2018: Copied from jrc3.m
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end
end %func