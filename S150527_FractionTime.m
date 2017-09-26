load D140330_Landmark;
load D140422_Controls;
%load img0_WP; %load probe trial background images

%% outlier fraction


%% Fig. 2C. Visit fraction
% lim = [0 30];   mode = 'time/visit';
% lim = [];   mode = 'time';
lim = [];   mode = 'visit';

rectCrop = [493 1083 312 902];
pSample = .9;
[RGB_1, mrGrid1] = mapSearchTime(subsample(vsTrialPool_E, pSample), [], mode, lim);
[RGB_2, mrGrid2] = mapSearchTime(subsample(vsTrialPool_L, pSample), [], mode, lim);
[RGB_3, mrGrid3] = mapSearchTime(subsample(vsTrialPool_P, pSample), [], mode, lim, img0_P);
% 
% 
% figure; suptitle([mode ', probe trials']); colormap gray;
% subplot 131; imshow(RGB_1); title('Early'); axis(rectCrop);
% subplot 132; imshow(RGB_2); title('Late');     axis(rectCrop);
% subplot 133; imshow(RGB_3); title('Probe');  axis(rectCrop);
% 

% Const
nGrid = 20;
thresh = 1/exp(1);
csPhase = {'Early', 'Late', 'Probe'};

% Input
csGrids = {mrGrid1, mrGrid2, mrGrid3};

% Range selection

cvVisit = cell(3,1);
for iPhase=1:numel(csPhase)
    xyCr = rectCrop / nGrid;
    viX = round(xyCr(1)):round(xyCr(2));
    viY = round(xyCr(3)):round(xyCr(4));

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

arrayfun(@(i)disp(zproptest2(cvVisit0{i}, cvVisit{i})),1:3)

%% fraction of time comparison
f1 = @(f, m1, m2) sum(f(m1)) / sum(f(m2));
figure;
mlF = getImageMask(img0, [0 15], 'FOOD');
mlA = getImageMask(img0, [0 50], 'CENTRE');
plotBarErr(f1(region_E, mlF, mlA), f1(region_L, mlF, mlA), f1(region_P, mlF, mlA));
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('Fraction of time near food (<15cm)');

%% time spent near landmarks or food (Fig. 

vsPhase = {'E', 'L', 'P'};
cvLM = cell(4,3);
cvLM1 = cell(3,1);

cviAnimal = cell(3,1);
cvrFracLm = cell(3,1);
cvrFracFood = cell(3,1);

for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    viAnimal = zeros(size(vsTrialPool));
    vrFracLm = zeros(size(vsTrialPool));
    vrFracFood = zeros(size(vsTrialPool));
    for iTrial = 1:numel(vsTrialPool)
        Strial = vsTrialPool(iTrial);
        viAnimal(iTrial) = Strial.iAnimal;
        S = poolTrials_location(Strial);
        vlLm = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
        vlFood = S.vrDf <= 15;
        vrFracLm(iTrial) = mean(vlLm(S.vlZone));
        vrFracFood(iTrial) = mean(vlFood(S.vlZone));
    end
    cviAnimal{iPhase} = viAnimal;
    cvrFracLm{iPhase} = vrFracLm;
    cvrFracFood{iPhase} = vrFracFood;
end
    


figure;
boxplot_cell(cvrFracLm);
title('Fraction of time near landmarks (<3cm)');
set(gca, 'YLim', [0 .4]);

figure;
boxplot_cell(cvrFracFood);
title('Fraction of time near Food (<15cm)');
set(gca, 'YLim', [0 .8]);


%%
[p, S] = zproptest2(vl1, vl2)

%% Fig. 3E. Search time fraction. kruskal-wallis test

kwtest_jjj(cvrFracLm);
ylim([0 .4]);
ylabel('Search time fraction');
%     1.0000    2.0000   18.9940   37.0139   55.0338    0.0000
%     1.0000    3.0000  -11.8622   12.8125   37.4872    0.5163
%     2.0000    3.0000  -45.3315  -24.2014   -3.0712    0.0186

kwtest_jjj(cvrFracFood);
ylim([0 .8]);
ylabel('Search time fraction');
%     1.0000    2.0000  -55.6244  -37.3472  -19.0701    0.0000
%     1.0000    3.0000  -67.2562  -42.2292  -17.2022    0.0002
%     2.0000    3.0000  -26.3138   -4.8819   16.5499    0.9293

