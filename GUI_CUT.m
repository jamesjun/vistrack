% input: iFrame, hObject

if isnan(iFrame) || isempty(iFrame), return; end
handles = guidata(hObject);

%Flip X, Y
%xy_names = {'CoM', 'head', 'head-mid', 'mid', 'tail-mid', 'tail'};
viF = iFrame+1:numel(handles.TC); %frame range to remove
if isempty(viF), return; end
handles.XC(viF, :) = [];
handles.YC(viF, :) = [];
handles.TC(viF) = [];
handles.XC_off(viF) = [];
handles.YC_off(viF) = [];

%ang_names = {'CoM', 'head-mid', 'tail-mid', 'body-bend', 'tail-bend'};
handles.AC(viF, :) = [];

handles.FLIM(2) = iFrame;
handles.MOV(:,:,viF) = [];

%update
guidata(hObject, handles);
msgbox(sprintf('Cut up to Frame %d', iFrame));

