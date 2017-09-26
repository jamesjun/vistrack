function [RGB, mrPlot] = mapSearchTime(vsTrial, strVar, mode, lim, img0)
nGrid = 20;
nTime = 25;
angRot = -1.1590; %deg

%---------------------
% Format input
if nargin < 2
    strVar = [];
end
if ~isempty(strVar)
    IPI = poolTrials_IPI(vsTrial);
    vrZ = getfield(IPI, strVar);
end
if nargin < 4
    lim = [];
end
if nargin < 5
    img0 = [];
end

%background image processing
xy0 = vsTrial(1).xy0;
if isempty(img0)
    img0 = vsTrial(1).img0;    
end
mlMask = getImageMask(img0, [0 60], 'CENTRE');
img0 = imrotate(imadjust(img0), angRot, 'nearest', 'crop');

%rotate vrX, vrY, and images
vrX = poolVecFromStruct(vsTrial, 'vrX');
vrY = poolVecFromStruct(vsTrial, 'vrY');
rotMat = rotz(angRot);    rotMat = rotMat(1:2, 1:2);
mrXY = [vrX(:) - xy0(1), vrY(:) - xy0(2)] * rotMat;
vrX = mrXY(:,1) + xy0(1);
vrY = mrXY(:,2) + xy0(2);

viX = ceil(vrX/nGrid);
viY = ceil(vrY/nGrid);
[h, w] = size(img0);
h = h / nGrid;
w = w / nGrid;
mnVisit = zeros(h, w);
mnTime = zeros(h, w);
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        viVisit = find(vlY & (viX == ix));        
        mnTime(iy,ix) = numel(viVisit);
        nRepeats = sum(diff(viVisit) < nTime); % remove repeated counts
        mnVisit(iy,ix) = numel(viVisit) - nRepeats;        
    end
end

mrTperV = mnTime ./ mnVisit;

if nargin < 3
    mode = 'time/visit';
end
switch lower(mode)
    case 'time'
        mrPlot = mnTime;
    case 'visit'
        mrPlot = mnVisit;
    case 'time/visit'
        mrPlot = mrTperV;
end


mnVisit1 = imresize(mrPlot, nGrid, 'nearest');
mnVisit1(~mlMask) = 0;

if isempty(lim)
    lim = [min(mnVisit1(:)) max(mnVisit1(:))];
end

mrVisit = uint8((mnVisit1 - lim(1)) / diff(lim) * 255);
RGB = rgbmix(img0, mrVisit, mlMask);

