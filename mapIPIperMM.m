function [mnVisit, mrVisit, mnI] = mapIPIperMM(S, img0, varPlot)
angRot = -1.1590; %deg
nGrid = 20;

vrX = S.vrX;
vrY = S.vrY;
if nargin < 3
    varPlot = 'vrD';
end
eval(sprintf('vrZ = S.%s;', varPlot)); 
if nargin < 2
    img0 = [];
end
if isempty(img0)
    img0 = S.img0;
end
img0 = imadjust(img0);
xy0 = S.xy0;
mlMask = getImageMask(img0, [0 60], 'CENTRE');

% rotational correction
img0 = imrotate(img0, angRot, 'nearest', 'crop');
rotMat = rotz(angRot);    rotMat = rotMat(1:2, 1:2);
mrXY = [vrX(:) - xy0(1), vrY(:) - xy0(2)] * rotMat;
vrX = mrXY(:,1) + xy0(1);
vrY = mrXY(:,2) + xy0(2);

% Grid counting
viX = ceil(vrX/nGrid);
viY = ceil(vrY/nGrid);
[h, w] = size(img0);
h = h / nGrid;
w = w / nGrid;
mnVisit = zeros(h, w);
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        mnVisit(iy,ix) = mean(vrZ(vlY & (viX == ix)));
    end
end
mnVisit = 1./mnVisit; %flip
mnVisit(isnan(mnVisit)) = 0;
mnVisit1 = imresize(mnVisit, nGrid, 'nearest');

mnVisit1(~mlMask) = 0;
lim = [0 1];
mrVisit = uint8((mnVisit1 - lim(1)) / diff(lim) * 255);
mnI = rgbmix(img0, mrVisit, mlMask);
