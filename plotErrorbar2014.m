function plotErrorbar(mrX, mlPath_E1, mlPath_L1, bootFcn, ystr)
%n by 3 (val, low, high)

[vrY, vrE, p] = bootCI_P(bootFcn, mrX(mlPath_E1), mrX(mlPath_L1));
% [bootFcn(mrX(mlPath_E1)), bootFcn(mrX(mlPath_L1))];

n = numel(vrY);
vrX = 1:n;
errorbar(vrX, vrY, vrE, 'r.'); hold on;
bar(vrX, vrY, .5);
set(gca, 'XLim', [.5 n+.5]);
    
set(gca, {'XTick', 'XTickLabel'}, {[1,2], {'Early', 'Late'}});
ylabel(ystr);

% [h p] = kstest2(mrX(mlPath_E1), mrX(mlPath_L1));
% [h p] = ttest2(mrX(mlPath_E1), mrX(mlPath_L1));

% title(sprintf('p = 10^{%f}', log10(p)));
if p < .0001
    title(sprintf('****, p = %f', p));
elseif p < .001
    title(sprintf('***, p = %f', p));
elseif p < .01
    title(sprintf('**, p = %f', p));
elseif p < .05
    title(sprintf('*, p = %f', p));
else
    title(sprintf('p = %f', p));
end

end


function [vrY, vrE, p] = bootCI_P(bootFcn, vrA, vrB)
nboot = 1000;
vrY = [bootFcn(vrA), bootFcn(vrB)];
[ciA, vrBootA] = bootci(nboot, {bootFcn, vrA});
[ciB, vrBootB] = bootci(nboot, {bootFcn, vrB});
vrE = abs(vrY - [ciA(1), ciB(1)]);
p=nan;
%p = ttest2_jjj(vrBootA, vrBootB, vrE(1), vrE(2), numel(vrA), numel(vrB));

% p = pval_kstest2(vrBootA, vrBootB, numel(vrA), numel(vrB));
% star rating
% * P ? 0.05 
% ** P ? 0.01 
% *** P ? 0.001 
% ****  P ? 0.0001 (see note) 


% http://faculty.psy.ohio-state.edu/myung/personal/course/826/bootstrap_hypo.pdf
% two-sample bootstrap hypothesis test
% tobs = vrY(1) - vrY(2);
% nboot1 = 3000;
% [bootstat, bootsam] = bootstrp(nboot1, bootFcn, [vrA(:); vrB(:)]);
% n = 0;
% for i=1:nboot1
%     t = bootFcn(bootsam(1:numel(vrA), i)) - bootFcn(bootsam(end-numel(vrB)+1:end, i));
%     n = n + (t > tobs);
% end
% p = (n / nboot1); %two-tailed

end %func