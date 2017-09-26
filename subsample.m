function [vr, vi] = subsample(vr, prob)
if prob>=1, return; end
if prob<=0, vr=[]; end

n = round(numel(vr) * prob);
% vi = randi(numel(vr), [n,1]);
vi = randperm(numel(vr), n);
% vi = round(linspace(1, numel(vr), n));
vr =vr(vi);