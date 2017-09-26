%% load data

load D140330_Landmark;
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csPlot = {'Early', 'Late', 'Probe'}; 

%% Fig. 3B.

%%

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
subplot 121;
boxplot_cell(cvrFracLm, 'mean-sem');
title('Fraction of time near landmarks (<3cm)');
set(gca, 'YLim', [0 .4]);

subplot 122;
boxplot_cell(cvrFracFood, 'mean-sem');
title('Fraction of time near Food (<15cm)');
set(gca, 'YLim', [0 .8]);


%% Fig. 4B. heading error angle

csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 

figure;
viZone = [1];
csZone = {'Az', 'Lm<3', 'F<15'};
for iZone1=1:numel(viZone)
    subplot(1,numel(viZone),iZone1);
    csZ = plotPhase(csTrials, viZone(iZone1), 'vrHeadingErr', []);
    xlabel('F<15, mean+-sem');
    ylabel('Heading error');
    title(csZone{iZone1});
    set(gca, 'YTick', 0:30:90);
    ylim([0 90]);
end

%%

csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csZ = plotPhase(csTrials, 1, 'vrHeadingErr', []);
figure; kwtest_jjj(csZ);
ylim([30 70]);
set(gca, 'YTick', 30:10:70);

ylabel('Heading error');
title('Error reduction with learning');
set(gca, 'XTick', 1:3, 'XTickLabel', {'Early', 'Late', 'Probe'});



%% Fig. 6B. Dist per DIPI
figure;
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csPlot = {'Early', 'Late', 'Probe'}; 

viZone = [1, 2, 3];
csZone = {'Az', 'Lm<3', 'F:3-15'};
for iZone1=1:numel(viZone)
    subplot(1,numel(viZone),iZone1);
    plotPhaseIpi(csTrials, viZone(iZone1), '1./nanmean(IPI.vrD(IPI.vlZ0))', []);
    %xlabel('F<15, mean+-sem');
    ylabel('Sampling density (cm^{-1})');
    title(csZone{iZone1});
    set(gca, 'YTick', 0:2:10);
    ylim([0 10]);
end

csZ = plotPhaseIpi(csTrials, viZone(2), '1./nanmean(IPI.vrD(IPI.vlZ0))', []);
figure; kwtest_jjj(csZ); title('landmark'); ylim([0 12]); set(gca, 'YTick', 0:3:12);

csZ = plotPhaseIpi(csTrials, viZone(3), '1./nanmean(IPI.vrD(IPI.vlZ0))', []);
figure; kwtest_jjj(csZ); title('near food'); ylim([0 12]); set(gca, 'YTick', 0:3:12);


%% backward swimming IPI
figure;
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csPlot = {'Early', 'Late', 'Probe'}; 

iZone = [3];
csTitle = {'Forward', 'Backward'};
csCmd = {'1./nanmean(IPI.vrD(IPI.vrV>0 & IPI.vlZ0))', '1./nanmean(IPI.vrD(IPI.vrV<0 & IPI.vlZ0))'};
for iCmd=1:numel(csCmd)
    subplot(1,numel(csCmd), iCmd);
    csZ = plotPhaseIpi(csTrials, iZone, csCmd{iCmd}, []);
    ylabel('Sampling density (cm^{-1})');
    title(csTitle{iCmd});
    set(gca, 'YTick', 0:10:40);
    ylim([0 40]);
end

%% prob. backward swimming

csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 

csZ = plotPhase(csTrials, 2, 'vrV<0', []);
figure; kwtest_jjj(csZ);
% ylim([30 70]);
% set(gca, 'YTick', 30:10:70);

ylabel('Prob. B-scans');
title('Near Landmark');
set(gca, 'XTick', 1:3, 'XTickLabel', {'Early', 'Late', 'Probe'});