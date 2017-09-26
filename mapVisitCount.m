function [RGB, mnVisit, mrTperV] = mapVisitCount(vsTrial, angRot, iAnimal)
nGrid = 20;
nTime = 25;

if nargin < 2
   angRot = [];
end
if nargin < 3
    iAnimal = [];
end

vrX = poolVecFromStruct(vsTrial, 'vrX', [], iAnimal);
vrY = poolVecFromStruct(vsTrial, 'vrY', [], iAnimal);
    
if ~isempty(angRot)
    angRot = -1.1590; %deg    
    xy0 = vsTrial(1).xy0;
    img0 = imadjust(vsTrial(1).img0);
    img0 = imrotate(img0, angRot, 'nearest', 'crop');
end

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
mnVisit1 = imresize(mnVisit, nGrid, 'nearest');
mrVisit = uint8(mnVisit1 / max(mnVisit1(:)) * 255);
mlMask = getImageMask(img0, [0 60], 'CENTRE');
RGB = rgbmix(img0, mnVisit1, mlMask);

% if nargout == 0
    
%     figure;
    imshow(RGB);
% end

mrTperV = mnTime ./ mnVisit;