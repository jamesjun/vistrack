function [mnVisit, mnVisit1] = calcGridStats(vsTrialPool, img0, varname, fun2, mlMask)
% pixpercm = 1053.28/(sqrt(2)*100);
% nGrid = 20; %2.6854cm/grid
nGrid = 25; %3.3567cm/grid

% nGrid = 25;
% nTime = 1; %20 msec
fEODAs = 0;

vrX = poolVecFromStruct(vsTrialPool, 'vrX');
vrY = poolVecFromStruct(vsTrialPool, 'vrY');
switch upper(varname)        
    case 'EODAS'
        vrZ = poolVecFromStruct(vsTrialPool, 'EODA');
        fEODAs = 1;
        fun1 = @(x)calcDistAsym(x);
    otherwise
        vrZ = poolVecFromStruct(vsTrialPool, varname);
        fun1 = @(x)mean(x);
end
if nargin >= 4
    fun1 = @(x)fun2(x);
end

%Averaging
% mnVisit = ...
%     calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, [0, 0]) * 1/3 + ...
%    (calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, [nGrid, 0]/2) + ...
%     calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, [-nGrid, 0]/2) + ...
%     calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, [0, nGrid]/2) + ...
%     calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, [0, -nGrid]/2))/6;

mnVisit = calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, [0, 0]);

mnVisit1 = imresize(mnVisit, nGrid, 'nearest');

% maxVal = nanstd(mnVisit1(~isinf(mnVisit1)))*2;
maxVal = 2;

disp(maxVal);
% maxVal = 10;

if nargout == 0
    mnVisit1(~mlMask) = 0;
    mrVisit = uint8(mnVisit1 / maxVal * 255);
    figure;
    
    if nargin >= 5
        imshow(rgbmix(img0, mrVisit, mlMask));
    else
        imshow(rgbmix(img0, mrVisit));
    end
end
end

function mnVisit = calcGrid(vrX, vrY, vrZ, fun1, img0, nGrid, xy0)

vrX = vrX + xy0(1);    
vrY = vrY + xy0(2);

viX = ceil(vrX/nGrid);
viY = ceil(vrY/nGrid);
[h, w] = size(img0);
h = h / nGrid;
w = w / nGrid;
mnVisit = zeros(h, w);
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        mnVisit(iy,ix) = fun1(vrZ(vlY & (viX == ix)));
    end
end
end