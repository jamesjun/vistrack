function [n, BW, BW1] = calcAreaFWHM(mrI)

[maxI, vi] = max(mrI(:));
[i,j] = ind2sub(size(mrI), vi);

BW = mrI >= maxI/exp(1);
BW1 = bwselect(BW, j, i, 8);
n = sum(BW1(:));

if nargout == 0
    figure; 
    subplot 121; 
    imagesc(BW);
    subplot 122; 
    imagesc(BW1);
end