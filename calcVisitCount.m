function [mnVisit1, mnVisit] = calcVisitCount(vsTrialPool, img0, mlMask, nGrid)
% pixpercm = 1053.28/(sqrt(2)*100);

if nargin < 4
    nGrid = 20; %3.3567cm/grid
end
% nGrid = 25; %2.6854cm/grid
nTime = 25; %250 msec

% nGrid = 25;
% nTime = 1; %20 msec

vrX = poolVecFromStruct(vsTrialPool, 'vrX');
vrY = poolVecFromStruct(vsTrialPool, 'vrY');

viX = ceil(vrX/nGrid);
viY = ceil(vrY/nGrid);
[h, w] = size(img0);
h = h / nGrid;
w = w / nGrid;

mnVisit = zeros(h, w);
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        viVisit = find(vlY & (viX == ix));        
        nRepeats = sum(diff(viVisit) < nTime); % remove repeated counts
        mnVisit(iy,ix) = numel(viVisit) - nRepeats;
    end
end

mnVisit1 = imresize(mnVisit, nGrid, 'nearest');

if nargin >= 3 && ~isempty(mlMask)
    mnVisit1(~mlMask) = 0;
    mnVisit(~imresize(mlMask, 1/nGrid)) = 0;
end

if nargout == 0
    mrVisit = uint8(mnVisit1 / max(mnVisit1(:)) * 255);
    
    figure;
%     imagesc(mrVisit);
    imshow(rgbmix(img0, mrVisit));
end
