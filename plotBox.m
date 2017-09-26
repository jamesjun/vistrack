function plotBox(vrX, vrY, vrEdges)
% [0,5,10,15,30]
vrQuantile = [0 .95 1]; %0:.2:1;
if nargin < 3
    vrEdges = quantile(vrX, vrQuantile); %not in a linear scale
end

labels = {};
% for i=1:(numel(vrEdges)-1)
%     labels{i} = sprintf('~%0.1f', vrEdges(i+1));
% end
for i=1:(numel(vrEdges)-1)
    labels{i} = sprintf('%0.0f%% (%0.2f)~', vrQuantile(i)*100, vrEdges(i));
end

% labels = {'<5', '5~10', '10~15', '>15'};
x_ordinal = ordinal(vrX, labels, [], vrEdges);

% figure; 
boxplot(vrY, x_ordinal);
L1 = x_ordinal == labels{1};
L2 = x_ordinal == labels{2};
[h,p] =kstest2(vrY(L1), vrY(L2));
title(sprintf('p=%f, n1=%d, n2=%d', p, sum(L1), sum(L2)));
% xlabel(str_x);
% ylabel(str_y);
% set(gca, 'YLim', [0 8], 'YTick', 0:4:8);
% grid on;