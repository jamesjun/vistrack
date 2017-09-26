function plotXcorr(X, Y)
subtMu = @(x)x - nanmean(x);
if nargin < 2
    Y = X;
end
n = 1000;

plotX = -n:n;
plotY = xcorr(subtMu(X), subtMu(Y), n, 'coeff');
plot(plotX, plotY,'r'); hold on;
plot(get(gca, 'XLim'), [0 0], 'k-');
plot(get(gca, 'XLim'), 1/exp(1)*[1 1], 'k-');
end