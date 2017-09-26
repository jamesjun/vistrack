function y = sem(vr)

y = std(vr) / sqrt(numel(vr));
