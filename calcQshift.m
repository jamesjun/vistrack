function vrQ = calcQshift(vrX)
% scales between -1 to 1 
n = numel(vrX);

[~, vi] = sort(vrX);
[~, vrQ] = sort(vi);
n0 = (n+1)/2; %mid value
n1 = (n-1)/2; %half value
vrQ = (vrQ - n0) / n1;
