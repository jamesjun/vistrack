function showMovieFrame(pos, handles, hImg, hFig, chLine, axImg)
persistent idxprev vrPlot fUsed;

if ~isempty(fUsed)
    return;
else
    fUsed = 1;
end
    
% figure(hFig);
tsel = pos(1);
idx = round(interp1(handles.TC, 1:numel(handles.TC), tsel));

fDraw = 1;
if ~isempty(idxprev)
    if idx == idxprev
        fDraw = 0;
    end
end
        
if fDraw
    title(axImg, sprintf('T_{ADC} = %0.3f; Frame# = %d', tsel, idx));
    set(hImg, 'CData', imadjust(handles.MOV(:,:,idx))); 
    for i=1:numel(vrPlot)
        try, delete(vrPlot(i)); catch, end;
    end
    % show marker
    Xc = handles.XC(idx,:);
    Yc = handles.YC(idx,:);
    xo = handles.XC_off(idx,:);
    yo = handles.YC_off(idx,:);
    hold(axImg, 'on');   
    vrPlot(1) = plot(axImg, Xc(2) - xo, Yc(2) - yo, 'ro');
    vrPlot(2) = plot(axImg, Xc(2:3) - xo, Yc(2:3) - yo, 'r-');
    vrPlot(3) = plot(axImg, Xc(5:6) - xo, Yc(5:6) - yo, 'r-');
end

for i=1:numel(chLine)
    pos1 = getPosition(chLine{i});
    pos1([1, 2]) = tsel;
    setPosition(chLine{i}, pos1);
end

idxprev = idx;
fUsed = [];
end