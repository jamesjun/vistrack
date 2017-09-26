function [vaf, Yf] = calcVAF(X, Y, vl)
%m: slope, b: intercept, corr
% VAF. variance accounted for

if nargin >=3
    X = X(vl);
    Y = Y(vl);
end
% fit1 = fit(X, Y,'poly1');
% Yf = feval(fit1, X);

if size(X,2) == 1 || size(X,1) == 1
    Yf = fitPoly1(X(:), Y(:));
else
    Yf = fitPoly1(X, Y);
%     Yf = feval(fit(X, Y,'poly11'), X);
end

vaf = 1 - nanvar(Yf - Y) / nanvar(Y);

