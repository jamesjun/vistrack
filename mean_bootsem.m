function mr = mean_bootsem(mrData)

nCol = size(mrData,2);
mr = zeros(nCol, 3);

for iCol = 1:nCol
    vrData = mrData(:,iCol);
    mr(iCol,1) = nanmean(vrData);
    mr(iCol,2:3) = bootci(1000, {@(y)nanmean(y), vrData});
end