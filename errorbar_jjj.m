function errorbar_jjj(vrX, mrIQR)

% if nargin<3, bin_width=.1; end
if isempty(vrX), vrX=1:size(mrIQR,1); end
errorbar(vrX, mrIQR(:,2), mrIQR(:,2)-mrIQR(:,1), mrIQR(:,3)-mrIQR(:,2), 'k-', 'MarkerSize', 1);
set(gca,'XLim', vrX([1,end]) + [-.5, .5]);