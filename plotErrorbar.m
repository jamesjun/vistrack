function plotErrorbar(mrX, mlPath_E1, mlPath_L1, bootFcn, ystr)
%n by 3 (val, low, high)

mr = [bootFcn(mrX(mlPath_E1)), bootFcn(mrX(mlPath_L1))];

n = size(mr, 2);
vrX = 1:n;
errorbar(vrX, mr(1,:), mr(1,:)-mr(2,:), mr(3,:)-mr(1,:), 'r.'); hold on;
bar(vrX, mr(1,:), .5);
set(gca, 'XLim', [.5 n+.5]);
    
set(gca, {'XTick', 'XTickLabel'}, {[1,2], {'Early', 'Late'}});
ylabel(ystr);

% [h p] = kstest2(mrX(mlPath_E1), mrX(mlPath_L1));
% title(sprintf('p = 10^{%f}', log10(p)));