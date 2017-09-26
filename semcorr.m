function semc = semcorr(vr)
% sem based on correlation time

vr = vr(~isnan(vr));
sd = std(vr);
n = numel(vr) / corrTime(vr);
semc = sd / sqrt(n);
end

function tau = corrTime(vrX)
thresh = 1/exp(1);
vrC = xcorr(vrX - mean(vrX), 'coeff');
vrC = vrC(ceil(end/2):end);
tau = find(vrC < thresh, 1, 'first');
end