function c = kwtest_jjj(cvr)

vr = [];
vi = [];
for i=1:numel(cvr)
    vr1 = cvr{i};
    vr = [vr; vr1(:)];
    vi1 = ones(size(vr1)) * i; %group
    vi = [vi; vi1(:)];
end

vlKill = isnan(vr);
vr(vlKill) = [];
vi(vlKill) = [];

[p, tbl, stats] = kruskalwallis(vr, vi, 'off');
boxplot(vr, vi);

c = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');
% c = multcompare(stats, 'Display', 'off');
% c = multcompare(stats);
disp(c);
