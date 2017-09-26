function [p, mx, my, semx, semy, nx, ny] = ttest2_jjj(x, y, bootFcn, ncrr)

if nargin < 3
    bootFcn = @mean;
end
if nargin < 4
    ncrr = [];
end

[mx, semx, nx] = calcBootSem(x, bootFcn, ncrr);
[my, semy, ny] = calcBootSem(y, bootFcn, ncrr);


tail = 0; %two sided

s2xbar = semx.^2;
s2ybar = semy.^2;
difference = mx - my;
dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
se = sqrt(s2xbar + s2ybar);
ratio = difference ./ se;


% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = 2 * tcdf(-abs(ratio),dfe);
elseif tail == 1 % right one-tailed test
    p = tcdf(-ratio,dfe);
elseif tail == -1 % left one-tailed test
    p = tcdf(ratio,dfe);
end

end %func