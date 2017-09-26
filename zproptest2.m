function [p, S] = zproptest2(vl1, vl2, n)
if nargin<3, n=2; end
k = nchoosek(n,2);
% based on http://en.wikipedia.org/wiki/Statistical_hypothesis_testing
% two-proportion z-test, pooled for H0: p1=p2

n1 = numel(vl1);
n2 = numel(vl2);
x1 = sum(vl1);
x2 = sum(vl2);
p = (x1+x2) / (n1+n2);

alpha = (x1/n1 - x2/n2) / sqrt(p * (1-p) * (1/n1 + 1/n2)); %alpha

p = normcdf(alpha);
if p > .5
    p = (1-p)*2;
else
    p = p * 2;
end
% p = 1-(1-p)^k;

[phat1, pci1] = binofit(x1, n1);
[phat2, pci2] = binofit(x2, n2);
S = struct('phat1', phat1, 'pci1', pci1, 'phat2', phat2, 'pci2', pci2);
%fprintf('p1: %f (%f-%f), p2: %f (%f-%f)\n', phat1, pci1, phat2, pci2);