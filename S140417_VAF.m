% VAF. variance accounted for

x = (1:10)';
y = 3*x+2;
y1 = y + rand(size(y)) * 10;

tic
fit1 = fit(x,y1,'poly1');
yf = feval(fit1, x);
toc
% 
% tic
% yf = fitPoly1(x, y1);
% toc


figure; plot(x,y1,'b', x,yf, 'r');

vaf = 1 - var(yf- y1) / var(yf)



calcVAF(x, y1)

%% 2D vaf
vaf = @(e,o)1 - var(e-o)/var(o);

x1 = rand(100,1);
x2 = rand(100,1);
y = 3*x1 + 2*x2 + 2;
y1 = y + rand(size(y)) * 1;

tic
vaf(feval(fit(x1, y1,'poly1'), x1), y1)
vaf(feval(fit(x2, y1,'poly1'), x2), y1)
vaf(feval(fit([x1,x2], y1,'poly11'), [x1,x2]), y1)
toc

calcVAF(x1, y1)
calcVAF(x2, y1)
calcVAF([x1,x2], y1)

tic
yf = fitPoly1(x, y1);
toc


calcVAF(x1, y1)
calcVAF(x2, y1)
calcVAF([x1(:),x2(:)], y1)