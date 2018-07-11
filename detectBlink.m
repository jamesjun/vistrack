function iFrame = detectBlink(handles, mode, fAsk)
% Returns the absolute frame
if nargin<3, fAsk = 1; end

mode = lower(mode);
vidobj = handles.vidobj;
switch mode
    case 'first'
    FLIM1 = [1, 300];
    FLIM1(2) = min(FLIM1(2), vidobj.NumberOfFrames);
    
    case 'last'
    FLIM1 = [-300, -1] + vidobj.NumberOfFrames + 1;
    FLIM1(1) = max(FLIM1(1), 1); 
    
    otherwise
    if isnumeric(mode)
        FLIM1 = [-75, 75] + mode;
        FLIM1 = min(max(FLIM1, 1), vidobj.NumberOfFrames); 
        mode = 'first';
    else
        error(sprintf('detectBlink-Unsupported mode-%s', mode));
    end
end

xyLED = handles.xyLED;

yRange = xyLED(2)+(-5:5);
xRange = xyLED(1)+(-5:5);
yRange1 = xyLED(2)+(-15:15);
xRange1 = xyLED(1)+(-15:15);

h=msgbox('Loading... (this will close automatically)');       
trImg = read(handles.vidobj, FLIM1);
try close(h); catch, end;

trImg = trImg(:,:,1,:);  

trInt = trImg(yRange, xRange,1,:);
trInt = mean(mean(trInt,1),2);
vrInt = trInt(:);
[vMax,iFrame] = max(vrInt);
[vMin,iFrameMin] = min(vrInt);
thresh = (vMax+vMin)/2; %vMax*.9 previously
iFrame = find(diff(vrInt>thresh)>0, 1, mode)+1;

if iFrame > 1 && iFrame < size(trImg,4)
    hfig = figure;     
    subplot 231; imshow(trImg(yRange1,xRange1,1,iFrame-1)); 
    title(num2str(iFrame-1));
    subplot 232; imshow(trImg(yRange1,xRange1,1,iFrame)); 
    title(num2str(iFrame));    
    subplot 233; imshow(trImg(yRange1,xRange1,1,iFrame+1)); 
    title(num2str(iFrame+1)); 
    subplot(2,3,4:6); bar(1:(diff(FLIM1)+1), vrInt); hold on;
    plot(iFrame*[1 1]+.5, get(gca, 'YLim'), 'r-');
    xlabel('Frame #'); ylabel('Intensity'); axis tight;
    set(gcf,'Name',handles.vidFname);
    
    button = questdlg('Is the blink detection correct?',handles.vidFname,'Yes','No','Yes');
    if strcmp(button, 'Yes')    
        fAskUser = 0;
    else
        fAskUser = 1;
    end
    try close(hfig), catch, end;
else
    fAskUser = 1;
end

if fAskUser && fAsk   
    implay(trImg);
    vcMsg = sprintf('Find the %s brightest blink, and close the movie', mode);
    uiwait(msgbox({vcMsg, handles.vidFname}));
    ans = inputdlg('Frame Number', 'First frame',1,{num2str(iFrame)});
    iFrame = str2num(ans{1});
    fprintf('Frame %d selected.\n', iFrame);
else
    iFrame = []; return
end

iFrame = iFrame + FLIM1(1) - 1;