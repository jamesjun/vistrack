function timerFcnPreview(hObject, event)

S=get(hObject, 'UserData');
%delete timer if figure closed
if ~ishandle(S.hFig)
    delete(hObject);
    return; 
end
handles = guidata(S.hObject);

S.iFrame = S.iFrame + S.REPLAY_STEP;
nFrames = min(size(handles.MOV,3), size(handles.XC,1));
if S.iFrame > nFrames
    if S.REPLAY_STEP == 1
        stop(hObject);
        return;
    else
        S.iFrame = nFrames;
    end
end

mrXC = bsxfun(@minus, handles.XC, handles.XC_off);
mrYC = bsxfun(@minus, handles.YC, handles.YC_off);
[~, vcDataID, ~] = fileparts(handles.vidFname);
set(S.hImg, 'CData', imadjust(handles.MOV(:,:,S.iFrame)));
XC = mrXC(S.iFrame,:);
YC = mrYC(S.iFrame,:);
nxy = numel(XC);
X1 = interp1(2:nxy, XC(2:end), 2:.1:nxy, 'spline');
Y1 = interp1(2:nxy, YC(2:end), 2:.1:nxy, 'spline');
set(S.hPlot(1), 'XData', XC(2), 'YData', YC(2));
set(S.hPlot(2), 'XData', XC(3:end), 'YData', YC(3:end));
set(S.hPlot(3), 'XData', X1, 'YData', Y1);
set(S.hPlot(4), 'XData', XC(1), 'YData', YC(1));
set(S.hTitle, 'String', ...
    sprintf('F1: %d; T1: %0.3f s, Step: %d (%s)', ...
        S.iFrame, S.TC1(S.iFrame), S.REPLAY_STEP, vcDataID));

set(hObject, 'UserData', S);