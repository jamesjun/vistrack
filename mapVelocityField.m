function [mnVisit, mrVisit, mnI] = mapVelocityField(S, img0)
angRot = -1.1590; %deg
nGrid = 20;

if nargin < 2
    img0 = [];
end

vrX = S.vrX;
vrY = S.vrY;
vrU = differentiate5(vrX);
vrV = differentiate5(vrY);

%direction correction
vl = vrU<0;
vrU(vl) = vrU(vl)*-1;
vrV(vl) = vrV(vl)*-1;

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
mrU = zeros(h, w);
mrV = zeros(h, w);
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        mrU(iy,ix) = mean(vrU(vlY & (viX == ix)));
        mrV(iy,ix) = mean(vrV(vlY & (viX == ix)));
    end
end
[mrX mrY] = meshgrid(1:h, 1:w);

if nargout == 0
    imshow(imadjust(img0)); hold on;
    quiver(mrX*nGrid - nGrid/2, mrY*nGrid - nGrid/2, mrU, mrV);
end