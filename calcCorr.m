function c = calcCorr(A,B, vl)
if nargin == 1
    B = A(:,2);
    A = A(:,1);    
end
if nargin >= 3
    A = A(vl);
    B = B(vl);
end

c = corrcoef(A,B);
c=c(2);