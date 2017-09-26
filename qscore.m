function vrQ = qscore(vrX)
% scales between -1 to 1 
n = numel(vrX);

[~, vi] = sort(vrX);
[~, vrQ] = sort(vi);
vrQ = (vrQ - 1) / (n-1);
