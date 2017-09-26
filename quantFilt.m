function vr = quantFilt(vr, quantLim)
qlim = quantile(vr(:), quantLim);
vr = vr(vr >= qlim(1) & vr < qlim(end));