function [maxX, areaX] = calcXcorr2Max(mrA, mrB, thresh)

if nargin < 3
    thresh = .6;
end
mrX = xcorr2(mrA, mrB);
maxX = max(max(mrX));
areaX = sum(sum(mrX > maxX * thresh));