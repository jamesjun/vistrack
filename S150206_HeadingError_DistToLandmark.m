% S140419 Heading direction analysis
load D140330_Landmark;


%% Heading angle grid map
% [mnVisit1, mnVisit] = calcGridStats(vsTrialPool_E, img0, 'EODAs');
CLIM = [0 90];
vsTrialPool = vsTrialPool_P;
titlestr = 'probe';
vcYlabel = 'survivor';

S=vsTrialPool(1);
img0 = S.img0;
xy0 = S.xy0;
calcD0 = @(x,y)sqrt((x-xy0(1)).^2 + (y-xy0(2)).^2) / pixpercm;        

%[vlZone, strZone] = getZone(S, iZone);

vlZone = logical(poolVecFromStruct(vsTrialPool, 'vlZone'));
vrAh = poolVecFromStruct(vsTrialPool, 'HANG');
vrX = poolVecFromStruct(vsTrialPool, 'vrX');
vrY = poolVecFromStruct(vsTrialPool, 'vrY');
xyf = [789, 681];
vrAf = cart2pol(xyf(1)-vrX, xyf(2)-vrY);
vrZ = abs(rad2deg(calcAngErr(vrAh, vrAf)));

figure; ksdensity(vrZ(vlZone), 'function', vcYlabel); 
xlabel('Heading error angle (deg)');
ylabel(vcYlabel);
title(titlestr);

%% distribution of heading error and compare E,L,P
vcYlabel = 'survivor'; %'survivor', 'pdf'
%vcZone = 'vlZone';  vcZoneDesc = 'active zone';
vcZone = 'vrD1<5';  vcZoneDesc = 'first landmark Anm2'; iAnimal = [];

titlestr = sprintf('Heading error angle distribution, near %s', vcZoneDesc);

S_P = poolTrials_RS(vsTrialPool_P, iAnimal);
S_E = poolTrials_RS(vsTrialPool_E, iAnimal);
S_L = poolTrials_RS(vsTrialPool_L, iAnimal);
figure; hold on;
vcPhase = 'ELP';
vcColor = 'gbr';
vhEntry = zeros(1,3);

for iPhase=1:3
    eval(sprintf('vrData = S_%c.vrHeadingErr(S_%c.%s);', ...
        vcPhase(iPhase), vcPhase(iPhase), vcZone));
    ecdf(vrData, 'function', 'survivor', 'bounds', 'on');
    vhLine = get(gca, 'children');
    for iLine=1:3
        set(vhLine(iLine), 'Color', vcColor(iPhase));
    end
    vhEntry(iPhase) = vhLine(3);
end

legend(vhEntry, {'Early', 'Late', 'Probe'}, 'box', 'off');
xlabel('Heading angle error (deg)');
ylabel(vcYlabel);
title(titlestr);
switch vcYlabel
    case 'survivor'
        set(gca, 'YTick', 0:.25:1); set(gca, 'YLim', [0 1]);
    case 'pdf'
        set(gca, 'YTick', 0:.005:0.02); set(gca, 'YLim', [0 .02]);   
end
set(gca, 'XLim', [0 90]);
set(gca, 'XTick', 0:15:90);


%% bar graph heading direction w.r.t the tangent of mid-point

warning off;
csCmd = {'abs(rad2deg(calcAngErr(RS.vrAh, RS.vrAf)))', '<|\theta_H - \theta_F|> (deg)', '@(x)median(x)', '@(x)sem(x)'; ...
         'abs(rad2deg(calcAngErr(RS.vrAh, RS.vrAt)))', '<|\theta_H - \theta_T|> (deg)', '@(x)median(x)', '@(x)sem(x)'; ...         
         'abs(rad2deg(calcAngErr(RS.vrAt, RS.vrAf)))', '<|\theta_T - \theta_F|> (deg)', '@(x)median(x)', '@(x)sem(x)'};     
   
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; csXlabel = {'Early', 'Late', 'Probe'};        
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csXlabel = {'Stable', 'None', 'Unstable'};
 
iZone = 3; %plot all active
plotAll(csTrials, csCmd, iZone, csXlabel);
set(gca, 'YLim', [0 90]);


%% pool backward swimming by trial
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 

figure;
viZone = [1, 2, 3];
csZone = {'Az', 'Lm<3', 'F<15'};
for iZone1=1:numel(viZone)
    subplot(1,numel(viZone),iZone1);
    plotPhase(csTrials, viZone(iZone1), 'vrHeadingErr', []);
    xlabel('F<15, mean+-sem');
    ylabel('Heading error');
    title(csZone{iZone1});
    set(gca, 'YTick', 0:30:90);
    ylim([0 90]);
end


%% bar graph plot of heading error and compare E,L,P
vcYlabel = 'survivor'; %'survivor', 'pdf'
%vcZone = 'vlZone';  vcZoneDesc = 'active zone';
vcZone = 'vrD1<5';  vcZoneDesc = 'second landmark Anm2'; iAnimal = 4;

titlestr = sprintf('Heading error angle distribution, near %s', vcZoneDesc);

S_P = poolTrials_RS(vsTrialPool_P, iAnimal);
S_E = poolTrials_RS(vsTrialPool_E, iAnimal);
S_L = poolTrials_RS(vsTrialPool_L, iAnimal);
figure; hold on;
vcPhase = 'ELP';
vcColor = 'gbr';
vhEntry = zeros(1,3);


%% example trajectory
vsTrialPool = vsTrialPool_P;
vcPhase = 'Probe';
img0 = imadjust(vsTrialPool(1).img0);

try
    figure(hFig); clf;
catch
    hFig = figure;
end

for iTrial=1:numel(vsTrialPool)
    Strial = vsTrialPool(iTrial);
    trial_PlotTraj(Strial, 'img0', img0);
    title(sprintf('%s trial: %d, %s', vcPhase, iTrial, Strial.dataID));
    pause(1);
end


%% 150208-1222: divide late trials to best and worst half and look at the landmark use
vsTrialPool = vsTrialPool_L;
iAnimal = 1;

vrPath = poolVecFromStruct(vsTrialPool, 'pathLen_cm', [], iAnimal);
vlTrialGood = vrPath < median(vrPath);

vsTrials_Good = vsTrialPool(vlTrialGood);
vsTrials_Bad = vsTrialPool(~vlTrialGood);

%S_good = poolTrials_RS(vsTrials_Good);
%S_bad = poolTrials_RS(vsTrials_Bad);

isNearLM = @(S,d) S.vrD1<d | S.vrD2<d | S.vrD3<d | S.vrD4<d;
%isNearLM = @(S,d) S.vrD1<d;
vrFracLm_good = zeros(size(vsTrials_Good));
for iTrial = 1:numel(vsTrials_Good)
    Strial = poolTrials_RS(vsTrials_Good(iTrial));
    vrFracLm_good(iTrial) = sum(isNearLM(Strial,3)) / sum(Strial.vlZone);
end
vrFracLm_bad = zeros(size(vsTrials_Bad));
for iTrial = 1:numel(vsTrials_Bad)
    Strial = poolTrials_RS(vsTrials_Bad(iTrial));
    vrFracLm_bad(iTrial) = sum(isNearLM(Strial,3)) / sum(Strial.vlZone);
end

figure; bar([mean(vrFracLm_good), mean(vrFracLm_bad)]); set(gca, 'XTickLabel', {'good', 'bad'});
ylabel('Fraction time spent near landmarks');
title(sprintf('Animal %d', iAnimal));


%% 150208-1503: example trajectories


