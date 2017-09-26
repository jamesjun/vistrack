function [mrPdf xb yb] = plotJointHist(Abin, Abin1, xb, yb)


dat = [Abin(:), Abin1(:)];
if nargin == 4
    n = hist3(dat,  {xb, yb});
else
    n = hist3(dat);
end
dx = mean(diff(xb));
dy = mean(diff(yb));

mrPdf = n'  / numel(Abin) / (dx * dy); %prob per area 
% n1( size(n,1) + 1 ,size(n,2) + 1 ) = 0; 

if nargin < 2
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
end

if nargout == 0
imagesc(mrPdf, 'xdata', xb, 'ydata', yb);
% imagesc(log(mrPdf)/log(10), 'xdata', xb, 'ydata', yb);
colormap hot; axis xy;
% axis equal;
end

% h = pcolor(xb,yb,log(n1)/log(10));
% h = pcolor(xb,yb,(n1));
% set(h, 'EdgeAlpha', 0);

