function vr = cell2vec(cm)
vr = [];
for i=1:numel(cm)
    vr1 = cm{i};
    vr = [vr; vr1(:)];
end