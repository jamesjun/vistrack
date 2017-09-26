% S140425_sampling density update

%% Dist per DIPI
% Dist per DIPI
% csTrials, csCmd, iZone, csX, nCols
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; csPlot = {'Early', 'Late', 'Probe'}; 
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csPlot = {'Stable', 'None', 'Unstable'}; 

csCmd = {'IPI.vrD', 'AZ', '@(x)1./mean(x)'; ...
         'IPI.vrD', 'LM', '@(x)1./mean(x)'; ...
         'IPI.vrD', 'NF', '@(x)1./mean(x)'; ...
         'IPI.vrD', 'F', '@(x)1./mean(x)'};
viZone = [1 2 3 4];

plotAll(csTrials, csCmd, viZone, csPlot);    

suptitle('Sampling density (cm^{-1})');


%% Dist per DIPI
% Dist per DIPI
% csTrials, csCmd, iZone, csX, nCols
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csX = {'getZone(IPI, 2)', 'Sampling density'; 'getZone(IPI, 3)', 'NF'; 'getZone(IPI, 4)', 'F'}; 
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csPlot = {'Stable', 'None', 'Unstable'}; 

csCmd = {'IPI.vrD', 'LM', '@(x)1./mean(x)'};
viZone = 1;

plotAll(csTrials, csCmd, viZone, csX);    

suptitle('Sampling density (cm^{-1})');

%% backward swimming

csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; csPlot = {'Early', 'Late', 'Probe'}; 
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csPlot = {'Stable', 'None', 'Unstable'}; 

csCmd = {'RS.vrV<0', 'AZ', '@(x)mean(x)'; ...
         'RS.vrV<0', 'LM', '@(x)mean(x)'; ...
         'RS.vrV<0', 'F<15', '@(x)mean(x)'};
viZone = [1 2 5];

plotAll(csTrials, csCmd, viZone, csPlot);    

suptitle('Prob. Backward swim');

%% backward swimming IPI

csTrials = {vsTrialPool_P}; 
csPlot = {'IPI.vrV>0', 'Forward'; 'IPI.vrV<0', 'Backward'}; 
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csPlot = {'Stable', 'None', 'Unstable'}; 

csCmd = {'IPI.vrD;', 'Sampling density', '@(x)1./mean(x)'};

viZone = [5];

plotAll(csTrials, csCmd, viZone, csPlot);    

suptitle('F<15:Probe');


%% DIPI analysis (E-Scan)

iPhase = 3;
iAnimal = 2;
figure;
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P};
csTitle = {'Early', 'Late', 'Probe'};
csColor = {'r', 'b', 'g'};

% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP};
% csTitle = {'Stable', 'None', 'Unstable'};

% iTrial = 2;
% vsTrial = csTrials{iPhase};
% S = vsTrial(iTrial);
% vrI = differentiate3(S.vtEOD);
% vrDI = differentiate3(vrI);

IPI = poolTrials_IPI(csTrials{iPhase}, iAnimal);
vrDI = IPI.vrDI * 1000;

vrX = (0:.000005:.0005)*1000;
vrDIp = vrDI(vrDI>=0);
vrDIn = -vrDI(vrDI<=0);
% vrDIp = hist(vrDIp,viX);
% vrDIn = hist(-vrDI(vrDI<=0),viX);
% figure; bar(viX, [vrDIp(:), vrDIn(:)]);
vpDIp = ksdensity(vrDIp, vrX, 'function', 'survivor');
vpDIn = ksdensity(vrDIn, vrX, 'function', 'survivor');
vpDIp(1)=1;
vpDIn(1)=1;

% figure; 
suptitle(csTitle{iPhase});
AX = [];

subplot 311; hold on; 
plot(vrX, vpDIp, csColor{iPhase}, vrX, vpDIn, csColor{iPhase});
legend('pos', 'neg');
ylabel('Survivor');
set(gca, 'YLim', [1e-4 1])
set(gca, 'YScale', 'log');
AX(end+1) = gca;

subplot 312; hold on;
plot(vrX, vpDIp./vpDIn, csColor{iPhase});
ylabel('neg/pos');
set(gca, 'YLim', [0 1.2])
set(gca, 'YTick', 0:.4:1.2)
AX(end+1) = gca;

subplot 313; hold on;
ylim = [min(vrSlope), max(vrSlope)];
vrSlope = differentiate5(vpDIp./vpDIn);
[~,imin] = min(vrSlope);
threshEscan = vrX(imin); %in msec
plot(vrX, vrSlope, csColor{iPhase});
set(gca, 'YLim', ylim);
plot(vrX(imin)*[1 1], ylim, 'k:');
ylabel('Slope');
% set(gca, 'YLim', [0 1.2])
% set(gca, 'YTick', 0:.4:1.2)
AX(end+1) = gca;

linkaxes(AX, 'x');
set(AX, 'XLim', [0 .2]);
set(AX, 'XTick', 0:.05:.2);
xlabel('|\DeltaIPI| (ms)'); 

% subplot 313; qqplot(vpDIp, vpDIn);
% set(gca, 'XLim', [0 .2])
% xlabel('\DeltaIPI (ms)'); ylabel('pos/neg');
% set(gca, 'XTick', 0:.05:.2);
% set(gca, 'YLim', [0 1.2])

%% E-Scan locations

IPI = poolTrials_IPI(csTrials{iPhase}, iAnimal);
vrDI = IPI.vrDI * 1000;
viTes = findDIsac(vrDI, -1*threshEscan); %E scan detection

% get distance distribution
vrX = ;
vrY = ;

