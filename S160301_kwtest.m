% S160301_KW test

load D140330_Landmark;
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 

nGrid = 20;
thresh = 1/exp(1);
csPhase = {'Early', 'Late', 'Probe'};
lim = [];   mode = 'visit';
rectCrop = [493 1083 312 902];

%%
xyCr = rectCrop / nGrid;
viX = round(xyCr(1)):round(xyCr(2));
viY = round(xyCr(3)):round(xyCr(4));

for iAnimal=1:4
    
    [~, cmrEarly{iAnimal}] = mapSearchTime(vsTrialPool_E, [], mode, lim, [], iAnimal);
    [~, cmrLate{iAnimal}] = mapSearchTime(vsTrialPool_L, [], mode, lim, [], iAnimal);
    [~, cmrProbe{iAnimal}] = mapSearchTime(vsTrialPool_P, [], mode, lim, [], iAnimal);
end

%%
% [~, mr] = mapSearchTime(vsTrialPool_E, [], mode, lim);
% [~, mr] = mapSearchTime(vsTrialPool_L, [], mode, lim);
[~, mr] = mapSearchTime(vsTrialPool_P, [], mode, lim);

mr = mr(viY, viX);
mean(mr(:) >= max(mr(:)) * thresh)
figure; imagesc(mr)

%% determine proportions

for iPhase = 1:3
    switch iPhase
        case 1, cmr = cmrEarly;
        case 2, cmr = cmrLate;
        case 3, cmr = cmrProbe;
    end
    for iAnimal=1:4
        mr0 = cmr{iAnimal};
        mrA = mr0(viY, viX); %active region cropped
        cutoff = max(mrA(:)) * thresh;
        fprintf('iAnimal:%d, %f\n', iAnimal, mean(mrA(:) >= cutoff));
    end
end

%% Const
    [~, mrGrid1] = mapSearchTime(vsTrialPool_E, [], mode, lim);
    [~, mrGrid2] = mapSearchTime(vsTrialPool_L, [], mode, lim);
    [~, mrGrid3] = mapSearchTime(vsTrialPool_P, [], mode, lim);

% Input
csGrids = {mrGrid1, mrGrid2, mrGrid3};

% Range selection

cvVisit = cell(3,1);
for iPhase=1:numel(csPhase)
    mr0 = csGrids{iPhase};
    mrA = mr0(viY, viX); %active region cropped
    cutoff = max(mrA(:)) * thresh;
    vl = mrA >= cutoff;
    cvVisit{iPhase} = double(vl(:)); 
end

figure; 
barplot_prop(cvVisit);
set(gca, 'XTickLabel', csPhase);
ylabel('Search area fraction');
% p(1,2) = 0.001115
% p(1,3) = 0.000020
% p(2,3) = 0.311572

zproptest2(cvVisit{1}, cvVisit{2})
zproptest2(cvVisit{2}, cvVisit{3})
zproptest2(cvVisit{1}, cvVisit{3})
