function [RGB, mnVisit] = mapTailBending(vsTrial, varStr)
nGrid = 20;
nTime = 25;
angRot = -1.1590; %deg

if nargin < 2
    varStr = 'vrDA'; %or varStr
end

vrX = poolVecFromStruct(vsTrial, 'vrX');
vrY = poolVecFromStruct(vsTrial, 'vrY');
vrZ = poolVecFromStruct(vsTrial, varStr);

xy0 = vsTrial(1).xy0;
img0 = imadjust(vsTrial(1).img0);
img0 = imrotate(img0, angRot, 'nearest', 'crop');

%rotate vrX, vrY, and images
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
% mnTime = zeros(h, w);
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        viVisit = (vlY & (viX == ix));        
%         mnTime(iy,ix) = numel(viVisit);
%         nRepeats = sum(diff(viVisit) < nTime); % remove repeated counts
%         mnVisit(iy,ix) = numel(viVisit) - nRepeats;
        mnVisit(iy,ix) = 1 ./ mean(abs(vrZ(viVisit) * 180 / pi));
    end
end

% Mix with RGB
mnVisit1 = imresize(mnVisit, nGrid, 'nearest');
mlMask = getImageMask(img0, [0 60], 'CENTRE');
mnVisit1(~mlMask) = 0;
lim = [0 2];
mrVisit = uint8((mnVisit1 - lim(1)) / diff(lim) * 255);
RGB = rgbmix(img0, mrVisit, mlMask);