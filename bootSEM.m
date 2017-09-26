function [mu, errL, errH] = bootSEM(vrY, fun1)

nBoot = 1000;
vec = bootci(nBoot, {@(y)mean(y), vrY}, 'type', 'cper');
kcorr = calcCorrTau(vrY); %correlation correction

mu = nanmean(vrY);

if nargin >= 2
    mu = fun1(mu);
    vec = sort(fun1(vec));
end

errL = (mu - vec(1)) * sqrt(kcorr);
errH = (vec(2) - mu) * sqrt(kcorr);