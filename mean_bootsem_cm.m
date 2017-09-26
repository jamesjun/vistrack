function [mr, cvData] = mean_bootsem_cm(cm)
nCol = size(cm,2);
mr = zeros(3,nCol);
cvData = cell(nCol,1);
for iCol = 1:nCol
    cvData{iCol} = cell2mat(cm(:,iCol));
    mr(:,iCol) = mean_bootsem(cvData{iCol});
end

mr = reshape(mr', [nCol,1,3]);

for i=1:nCol-1
    for j=i+1:nCol
%         kstest2_disp(cvData{i}, cvData{j}, sprintf('%d-%d', i,j));
    end
end

kwtest_jjj(cvData);