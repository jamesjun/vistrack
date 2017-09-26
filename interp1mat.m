function mrYi = interp1mat(vrX, mrY, vrXi, method)

if nargin < 4
    method = 'spline';    
end

if numel(vrX) ~= size(mrY,1)
    fFlip = 1;
else
    fFlip = 0;
end

if fFlip
    mrY = mrY';    
end
n = numel(vrXi);
m = size(mrY, 2);
mrYi = zeros(n, m);

for i=1:m
    mrYi(:,i) = interp1(vrX, mrY(:,i), vrXi, method);
end

if fFlip
    mrYi = mrYi';    
end
