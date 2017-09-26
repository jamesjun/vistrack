%% Grid based amp analysis

load D140330_Landmark;
load D140422_Controls;
load img0_WP; %load probe trial background images

% vsTrialPool_NP vsTrialPool_RP vsTrialPool_WP;
 
%% Fig. 10A, B Search time per grid

lim = [0 30];   mode = 'time/visit';
% lim = [];   mode = 'time';
% lim = [];   mode = 'visit';

rectCrop_P = [493 1083 312 902];
dxy = vsTrialPool_NP(1).xy0 - vsTrialPool_P(1).xy0;
rectCrop_NP = rectCrop + [dxy(1), dxy(1), dxy(2), dxy(2)];
dxy = vsTrialPool_WP(1).xy0 - vsTrialPool_P(1).xy0;
rectCrop_WP = rectCrop + [dxy(1), dxy(1), dxy(2), dxy(2)];

[RGB_1, mrGrid1] = mapSearchTime(vsTrialPool_P, [], mode, lim);
[RGB_2, mrGrid2] = mapSearchTime(vsTrialPool_NP, [], mode, lim);
[RGB_3, mrGrid3] = mapSearchTime(vsTrialPool_WP, [], mode, lim, img0_WP_A);

figure; suptitle([mode ', probe trials']); 
subplot 131; imshow(RGB_1); title('Landmark'); axis(rectCrop_P);
subplot 132; imshow(RGB_2); title('Empty');     axis(rectCrop_NP);
subplot 133; imshow(RGB_3); title('Random');  axis(rectCrop_WP);


%% Fig. 10B. Fraction of 1/e

% Const
nGrid = 20;
thresh = 1/exp(1);
csPhase = {'Stable', 'None', 'Unstable'};

% Input
csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP};
csGrids = {mrGrid1, mrGrid2, mrGrid3};

% Range selection
xyCr = rectCrop / nGrid;
viX = round(xyCr(1)):round(xyCr(2));
viY = round(xyCr(3)):round(xyCr(4));

cvVisit = cell(3,1);
for iPhase=1:numel(csTrials)
    mr0 = csGrids{iPhase};
    mrA = mr0(viY, viX); %active region cropped
    cutoff = max(mrA(:)) * thresh;
    vl = mrA >= cutoff;
    cvVisit{iPhase} = double(vl(:)); 
end

figure; plotCellErrorbar(cvVisit, [], [], 1);
set(gca, 'XTickLabel', csPhase);
ylabel('Fraction of search area');
axis([.5 3.5 0 1]); 
set(gca, 'YTick', 0:.2:1);
title('Search area (>1/e) as a fraction of active zone');


%% Fig. 10A, IPI/mm map <mm/IPI>^-1 per grid
rectCrop_P = [493 1083 312 902];

titleStr = 'IPI/mm, Late trials';

csPhase = {'Stable', 'None', 'Unstable'};
csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP};
img0 = [];
figure; 
for iPhase = 1:3
    subplot(1,3,iPhase);
    vsTrial = csTrials{iPhase};
%     S = poolDIPIperDist(vsTrial, []); %replaced by poolTrialsIPI
    S = poolTrials_IPI(vsTrial);
    switch iPhase
        case 1, rectCrop=rectCrop_P; img0=[];
        case 2, rectCrop=rectCrop_NP; img0=[];
        case 3, rectCrop=rectCrop_WP; img0=img0_WP_A;
    end
    [mnVisit, mnVisit1, mnI] = mapIPIperMM(S, img0);
    imshow(mnI);
    axis(rectCrop);
    title(csPhase(iPhase));
end
suptitle('<mm/IPI>^-1');


%% back-swim prob. box plot. distribution box plot
% strVar = 'RS.vpBack'; fun1 = [];
% csCmd = 'RS.vrV'; fun1 = @(x)abs(x);
csCmd = 'IPI.vrPmm'; fun1 = [];
iAnimal = 4;
csPlot = cell(3,1);
vrPlot = [];
vrGroup = [];
vrMu = zeros(3,1);
vrSem = zeros(3,1);
for iPhase = 1:3
    vsTrial = csTrials{iPhase};
    RS = poolTrials_RS(vsTrial, iAnimal);
    IPI = poolTrials_IPI(vsTrial, iAnimal);
    
    eval(sprintf('vrZ = %s;', csCmd));
    if ~isempty(fun1)
        vrZ = fun1(vrZ);
    end
    vrMu(iPhase) = mean(vrZ);
    vrSem(iPhase) = sem(vrZ);
%     csPlot{iPhase} = S.vpBack;
%     vrPlot = [vrPlot; S.vpBack];
%     vrGroup = [vrGroup; iPhase * ones(size(S.vpBack))];
end
% figure; boxplot(vrPlot, vrGroup, 'width', .75);
figure; plotBarError(vrMu, vrSem, [], {'Stable', 'None', 'Unstable'});
ylabel(csCmd);


% set(gca, 'YTick', 0:.4:1.2)

