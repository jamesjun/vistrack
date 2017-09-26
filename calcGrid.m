function [mnVisit, mnVisit1] = calcGrid(vrX, vrY, vrZ, img0, nGrid, mlMask)

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

mnVisit1 = imresize(mnVisit, nGrid, 'nearest');

if nargout == 0    
    disp(maxVal);
    mnVisit1(~mlMask) = 0;
    maxVal = quantile(mnVisit1(:), .95);
    mrVisit = uint8(mnVisit1 / maxVal * 255);
    figure;
    if nargin >= 6
        imshow(rgbmix(img0, mrVisit, mlMask));
    else
        imshow(rgbmix(img0, mrVisit));
    end
end