function boxplot_cell(cvr, vcMode)

if nargin < 2
    vcMode = 'median-iqr'; %or 'mean-sem'
end

nGroups = numel(cvr);
vrData = [];
mrData = zeros(3,nGroups);
for iGroup = 1:nGroups
    n = numel(cvr{iGroup});
    vrData = toVec(cvr{iGroup});
    switch lower(vcMode)
        case 'mean-sem'
            vec = bootci(1000, {@(y)nanmean(y), vrData});
            mrData(:,iGroup) = [nanmean(vrData), vec'];
        case 'median-iqr'
            mrData(:,iGroup) = quantile(vrData, [.5 .25 .75]);
    end
    hold on;
    %plot(iGroup, vrData, 'k.');
end

errorbar(1:nGroups, mrData(1,:), mrData(1,:)-mrData(2,:), mrData(3,:)-mrData(1,:),'ro');

set(gca, 'XTick', 1:nGroups);

%boxplot(vrData, viGroup,  'medianstyle', 'target', ...
%    'notch', 'off', 'whisker', 0, 'boxstyle', 'filled', 'symbol', '');

for i=1:(numel(cvr)-1)
    vl1 = cvr{i};
    for j=i+1:numel(cvr)
        vl2 = cvr{j};
        [h, p] = kstest2(vl1(:), vl2(:));
        fprintf('p(%d,%d) = %f, n1=%d, n2=%d\n', i, j, p, numel(vl1), numel(vl2));
    end
end
