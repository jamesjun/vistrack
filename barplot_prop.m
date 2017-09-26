function barplot_prop(cvl)

vrProp = zeros(numel(cvl), 3);

for i=1:numel(cvl)
    [phat, pci] = binofit(sum(cvl{i}), numel(cvl{i}));
    vrProp(i,:) = [phat, pci];
end

hold on;
vrXp = 1:numel(cvl);
bar(vrXp, vrProp(:,1));
errorbar(vrXp, vrProp(:,1), vrProp(:,3)-vrProp(:,1), vrProp(:,1)-vrProp(:,2), '.');

for i=1:numel(cvl)
    vl1 = cvl{i};
    for j=i+1:numel(cvl)
        vl2 = cvl{j};
        [p, S] = zproptest2(vl1, vl2);
        fprintf('p(%d,%d) = %f\n', i, j, p);
    end
end
set(gca, 'XTick', vrXp);