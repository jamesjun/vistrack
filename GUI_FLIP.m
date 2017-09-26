% input: iFrame, hObject

if isnan(iFrame) || isempty(iFrame), return; end
handles = guidata(hObject);

%Flip X, Y
%xy_names = {'CoM', 'head', 'head-mid', 'mid', 'tail-mid', 'tail'};
viF = iFrame:numel(handles.TC);
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
