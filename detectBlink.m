function [iFrame, xyLED] = detectBlink(handles, mode, fAsk)
% Returns the absolute frame
if nargin<3, fAsk = 1; end

mode = lower(mode);
vidobj = handles.vidobj;
switch mode
    case 'first'
    FLIM1 = [1, 300];
    FLIM1(2) = min(FLIM1(2), vidobj.NumberOfFrames);
    xyLED = []; % auto-detect
    
    case 'last'
    FLIM1 = [-300, -1] + vidobj.NumberOfFrames + 1;
    FLIM1(1) = max(FLIM1(1), 1); 
    xyLED = handles.xyLed;
    
    otherwise
    if isnumeric(mode)
        FLIM1 = [-75, 75] + mode;
        FLIM1 = min(max(FLIM1, 1), vidobj.NumberOfFrames); 
        mode = 'first';
    else
        error(sprintf('detectBlink-Unsupported mode-%s', mode));
    end
end

h=msgbox('Loading... (this will close automatically)', 'detect LED blink'); drawnow;    
% trImg = read(handles.vidobj, FLIM1);
% trImg = squeeze(trImg(:,:,1,:));
trImg = vid_read(handles.vidobj, FLIM1(1):FLIM1(2));
try close(h); catch, end;

if isempty(xyLED), xyLED = find_mov_max_(trImg); end
yRange = xyLED(2)+(-5:5);
xRange = xyLED(1)+(-5:5);
yRange1 = xyLED(2)+(-15:15);
xRange1 = xyLED(1)+(-15:15);


trInt = trImg(yRange, xRange,:);
trInt = mean(mean(trInt,1),2);
vrInt = trInt(:);
[vMax,iFrame] = max(vrInt);
[vMin,iFrameMin] = min(vrInt);
thresh = (vMax+vMin)/2; %vMax*.9 previously
iFrame = find(diff(vrInt>thresh)>0, 1, mode)+1;

if iFrame > 1 && iFrame < size(trImg,3)
    hfig = figure;     
    subplot 231; imshow(trImg(yRange1,xRange1,iFrame-1)); 
    title(num2str(iFrame-1));
    subplot 232; imshow(trImg(yRange1,xRange1,iFrame)); 
    title(num2str(iFrame));    
    subplot 233; imshow(trImg(yRange1,xRange1,iFrame+1)); 
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
        iFrame = nan;
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
end

if ~isnan(iFrame)
    iFrame = iFrame + FLIM1(1) - 1;
end
end %func


%--------------------------------------------------------------------------
function xyLed = find_mov_max_(trImg)
img_pp = (max(trImg,[],3) - min(trImg,[],3));
[~,imax_pp] = max(img_pp(:));
[yLed, xLed] = ind2sub(size(img_pp), imax_pp);
xyLed = [xLed, yLed];
end %func