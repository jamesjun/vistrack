function kstest2_disp(vr1, vr2, vc)

[h, p] = kstest2(vr1, vr2);
fprintf('p(%s)=%f, n1=%d, n2=%d\n', vc, p, numel(vr1), numel(vr2));