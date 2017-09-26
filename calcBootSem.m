function [mux, semx, nx, ncrr] = calcBootSem(x, bootFcn, ncrr)
if nargin < 3
    ncrr = [];
end
if isempty(ncrr)
    ncrr = calcCorrTau(x);
end
nx = numel(x) / ncrr;
nboot = 1000;
cix = bootci(nboot, {bootFcn, x}, 'type', 'norm'); %subsample
semx = diff(cix) / 2 * sqrt(ncrr);
mux = bootFcn(x);
end
