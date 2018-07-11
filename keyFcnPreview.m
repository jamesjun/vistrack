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
        csHelp = {'SPACE: start and stop', 'LEFT/RIGHT: backward/forward', ...
            'HOME/END: start/end', 'F: Flip head/tail', ...
            'UP/DOWN: Speed up/down', 'C: Trim video'};
        msgbox(csHelp);
        
    case 'space'  % toggle start and stop  
        if fRunning, stop(timer1);
        else start(timer1); end
    
    case {'leftarrow', 'rightarrow', 'f', 'home', 'end'}
    %'f' for flip, 'left' for back, 'right' for forward
        if fRunning, stop(timer1); end
        switch event.Key
            case 'leftarrow', S.iFrame = max(1, S.iFrame-S.REPLAY_STEP);
            case 'rightarrow', S.iFrame = min(nFrames, S.iFrame+S.REPLAY_STEP);
            case 'home', S.iFrame = 1;
            case 'end', S.iFrame = nFrames;
            case 'f', GUI_FLIP; %flip orientation
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