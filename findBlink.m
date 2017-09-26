function [viBlink, vrBlinkInt] = findBlink(vidobj, SYNC_PERIOD, iF0, FPS)

xyLED = [1555, 615];    %LED location
viRange = [-15 15];   %Temporal range
threshRatio = .9;
% FPS = vidobj.NumberOfFrames / vidobj.Duration;

yRange = xyLED(2)+(-5:5);
xRange = xyLED(1)+(-5:5);

viBlink = [];
vrBlinkInt = [];
iF = iF0;
h=msgbox('Sync in progress... (this will close automatically)');       
try
    while 1
        FLIM1 = iF + viRange;
        trImg = read(vidobj, FLIM1);
        trInt = trImg(yRange, xRange,1,:);
        trInt = mean(mean(trInt,1),2);
        vrInt = trInt(:);
        [vmax, imax] = max(vrInt);
        vmin = min(vrInt);
        ithresh = find(vrInt > vmax *threshRatio, 1, 'first');
        viBlink(end+1) = imax + FLIM1(1) - 1; %or use ithresh crossing
        vrBlinkInt(end+1) = vrInt(imax)/vmin;

        %next
%         if numel(viBlink) > 1
%             FPS = (viBlink(end)-viBlink(end-1)) / SYNC_PERIOD;
%         end
        iF = iF + round(SYNC_PERIOD * FPS);        
    end
catch
    ;%disp(lasterr);
end

try close(h); catch, end;