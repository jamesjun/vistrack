%S140417 variance accounted for
load D140330_Landmark;
load D140422_Controls;


%%
sem = @(x)std(x)/numel(x);
% csCmd = {'abs(IPI.vrD)', '<\Deltad>^{-1} (N/mm)'; ...
%          'abs(rad2deg(IPI.vrDHA))', '<\Delta\theta_H>^{-1} (N/deg)'; ...
%          'abs(rad2deg(IPI.vrDA))', '<\Delta\theta_T>^{-1} (N/deg)'; ...
%          'abs(RS.vrV)', 'Speed (cm/s)';...
%          'abs(rad2deg(RS.vrHAV))', '<\omega_H> (deg/s)';...
%          'abs(rad2deg(RS.vrV))', '<\theta_T> (deg/s)'};

csCmd = {'abs(IPI.vrD*10)', '<\Deltad> (mm/IPI)'; ...
         'abs(RS.vrV)', '<Speed> (cm/s)';...
         'abs(rad2deg(RS.vrAt))', '<\theta_T> (deg/s)'; ...
         'double(RS.vrV<0)', 'Prob. Back-swim'};

vcShape = 'o^sd'; %animal's shape
hfig = figure; 
nPlot = size(csCmd,1);     
for iCmd = 1:nPlot
         
[AX, AX1, csZ_V] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,1}, [], csCmd{iCmd,2});
close(gcf);
figure(hfig);    
subplot(1,nPlot,iCmd);
hold on;

for iAnimal = 1:4
    vrX = [1,2,3];
    vrY2 = [mean(csZ_V{iAnimal+1, 1, 1}), mean(csZ_V{iAnimal+1, 1, 2}), mean(csZ_V{iAnimal+1, 1, 3})];
    vrE2 = [sem(csZ_V{iAnimal+1, 1, 1}), sem(csZ_V{iAnimal+1, 1, 2}), sem(csZ_V{iAnimal+1, 1, 3})];
    errorbar(vrX, vrY2, vrE2, ['k' vcShape(iAnimal), ':']);
    vrY3 = [mean(csZ_V{iAnimal+1, 2, 1}), mean(csZ_V{iAnimal+1, 2, 2}), mean(csZ_V{iAnimal+1, 2, 3})];
    vrE3 = [sem(csZ_V{iAnimal+1, 2, 1}), sem(csZ_V{iAnimal+1, 2, 2}), sem(csZ_V{iAnimal+1, 2, 3})];
    errorbar(vrX, vrY3, vrE3, ['r' vcShape(iAnimal), ':']);
end
set(gca, 'XTick', 1:3);
set(gca, 'XLim', [.5 3.5]);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});

ylabel(csCmd{iCmd,2});

end

%% Corr stats 
csCmd = {'IPI.vrCorrV_D', 'IPI.vrCorrI_D', 'Corr. w/ \Deltad'; ...
         'IPI.vrCorrV_I', 'IPI.vrCorrA_I', 'Corr. w/ IPI';...
         'IPI.vrCorrVHA_A', 'IPI.vrCorrD_A', 'Corr. w/ \theta_T'};
sem = @(x)std(x)/numel(x);

vcShape = 'o^sd'; %animal's shape
hfig = figure; 
     
for iCmd = 1:size(csCmd,1)
         
[AX, AX1, csZ_V1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,1}, [], csCmd{iCmd,3});  close(gcf);
[AX, AX1, csZ_V2] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,2}, [], csCmd{iCmd,3});  close(gcf);

figure(hfig); subplot(1,3,iCmd); hold on;
for iAnimal = 1:4
    vrX = [1,2,3];
    vrY1 = [mean(csZ_V1{iAnimal+1, 1, 1}), mean(csZ_V1{iAnimal+1, 1, 2}), mean(csZ_V1{iAnimal+1, 2, 3})];
    vrE1 = [sem(csZ_V1{iAnimal+1, 1, 1}), sem(csZ_V1{iAnimal+1, 1, 2}), sem(csZ_V1{iAnimal+1, 2, 3})];
    errorbar(vrX, vrY1, vrE1, ['k' vcShape(iAnimal), ':']);
    vrY2 = [mean(csZ_V2{iAnimal+1, 1, 1}), mean(csZ_V2{iAnimal+1, 1, 2}), mean(csZ_V2{iAnimal+1, 3, 3})];
    vrE2 = [sem(csZ_V2{iAnimal+1, 1, 1}), sem(csZ_V2{iAnimal+1, 1, 2}), sem(csZ_V2{iAnimal+1, 3, 3})];
    errorbar(vrX, vrY2, vrE2, ['b' vcShape(iAnimal), ':']);
end
set(gca, 'XTick', 1:3);
set(gca, 'XLim', [.5 3.5]);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});

ylabel(csCmd{iCmd,3});

end


%% Vaf stats 
csCmd = {'IPI.vrVafV_D', 'IPI.vrVafI_D', 'Vaf w/ \Deltad'; ...
         'IPI.vrVafV_I', 'IPI.vrVafA_I', 'Vaf w/ IPI';...
         'IPI.vrVafVHA_A', 'IPI.vrVafV_A', 'Vaf w/ \theta_T'};
sem = @(x)std(x)/numel(x);

vcShape = 'o^sd'; %animal's shape
hfig = figure; 
     
for iCmd = 1:size(csCmd,1)
         
[AX, AX1, csZ_V1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,1}, [], csCmd{iCmd,3});  close(gcf);
[AX, AX1, csZ_V2] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,2}, [], csCmd{iCmd,3});  close(gcf);

figure(hfig); subplot(1,3,iCmd); hold on;
for iAnimal = 1:4
    vrX = [1,2,3];
    vrY1 = [mean(csZ_V1{iAnimal+1, 1, 1}), mean(csZ_V1{iAnimal+1, 1, 2}), mean(csZ_V1{iAnimal+1, 1, 3})];
    vrE1 = [sem(csZ_V1{iAnimal+1, 1, 1}), sem(csZ_V1{iAnimal+1, 1, 2}), sem(csZ_V1{iAnimal+1, 1, 3})];
    errorbar(vrX, vrY1, vrE1, ['k' vcShape(iAnimal), ':']);
    vrY2 = [mean(csZ_V2{iAnimal+1, 1, 1}), mean(csZ_V2{iAnimal+1, 1, 2}), mean(csZ_V2{iAnimal+1, 1, 3})];
    vrE2 = [sem(csZ_V2{iAnimal+1, 1, 1}), sem(csZ_V2{iAnimal+1, 1, 2}), sem(csZ_V2{iAnimal+1, 1, 3})];
    errorbar(vrX, vrY2, vrE2, ['b' vcShape(iAnimal), ':']);
end
set(gca, 'XTick', 1:3);
set(gca, 'XLim', [.5 3.5]);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});

ylabel(csCmd{iCmd,3});

end


%% Backward-swim stats. forward vs. backward
% 'double(RS.vrV<0)', 'Prob. Backswim'

csCmd = {'IPI.vrD*10; vlZ=vlZ&IPI.vrV>0;', 'IPI.vrD*10; vlZ=vlZ&IPI.vrV<0;', '<\Deltad> (mm/IPI)'; ...
         'abs(rad2deg(RS.vrHAV)); vlZ=vlZ&RS.vrV>0;', 'abs(rad2deg(RS.vrHAV)); vlZ=vlZ&RS.vrV<0;', '<\omega_H> (deg/s)'; ...
         'abs(RS.vrV); vlZ=vlZ&RS.vrV>0;', 'abs(RS.vrV); vlZ=vlZ&RS.vrV<0;', '<Speed> (cm/s)'; ...
         'IPI.vrI*1000; vlZ=vlZ&IPI.vrV>0;', 'IPI.vrI*1000; vlZ=vlZ&IPI.vrV<0;', '<IPI> (ms)'};

sem = @(x)std(x)/numel(x);

vcShape = 'o^sd'; %animal's shape
vcPhase = 'rbg';
hfig = figure; 
iZone = 1;     
vsTrialPool = [vsTrialPool_E, vsTrialPool_L, vsTrialPool_P];
nPlot = size(csCmd,1);
for iCmd = 1:nPlot
   
[AX, AX1, csZ_V1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,1}, [], csCmd{iCmd,3});  close(gcf);
[AX, AX1, csZ_V2] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    csCmd{iCmd,2}, [], csCmd{iCmd,3});  close(gcf);

figure(hfig); subplot(1,nPlot,iCmd); hold on;
vrX = [1,2];
for iPhase = 1:3
    for iAnimal = 1:4        
        vrY = [mean(csZ_V1{iAnimal+1, iZone, iPhase}), mean(csZ_V2{iAnimal+1, iZone, iPhase})];
        vrE = [sem(csZ_V1{iAnimal+1, iZone, iPhase}), sem(csZ_V2{iAnimal+1, iZone, iPhase})];
        errorbar(vrX, vrY, vrE, [vcPhase(iPhase) vcShape(iAnimal), ':']);
    end
end
set(gca, 'XTick', 1:2);
set(gca, 'XLim', [.5 2.5]);
set(gca, 'XTickLabel', {'Forward', 'Back'});

ylabel(csCmd{iCmd,3});

end


%% Approach vs. departure
% 'double(RS.vrV<0)', 'Prob. Backswim'

% csCmd = {'IPI.vrD*10', '<\Deltad> (mm/IPI)'; ... %R1
%          'abs(IPI.vrV)', '<Speed> (cm/s)'; ...
%          'abs(rad2deg(IPI.vrAh-IPI.vrAt))', '<\Deltatheta_E> (deg/s)'; ...
%          'IPI.vrIz', '<IPIz>'; ... %R2
%          'IPI.vrIq', '<IPIq>'; ...
%          'IPI.vrI*1000', '<IPI> (ms)'; ...
%          'IPI.vrWt', '<\omega_T> (ms)'; ... %R3
%          'IPI.vrWh', '<\omega_H> (ms)'; ...
%          'IPI.vrWv', '<\omega_V> (ms)'; ...
%          'IPI.vrAt', '<\theta_T> (ms)'};
     
csCmd = {'IPI.vrD*10', '<\Deltad> (mm/IPI)'; ... %R1        
         'abs(IPI.vrV)', '<Speed> (cm/s)'; ...         
         'IPI.vrAcc', '<Acc> (cm/s^2)'; ...
         'IPI.vrI*1000', '<IPI> (ms)'};     
iZone = 1; %plot all active
plotAll({vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}, csCmd, iZone, {'vlApp(IPI)', 'App'; 'vlDep(IPI)', 'Dep'});

%% search-time near food
csCmd = {'getZone(RS,2)', 'Prob. near Landmarks'; ... %R1        
         'getZone(RS,5)', 'Prob. near Food'};     
iZone = 1; %plot all active
plotAll({vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}, csCmd, iZone);

% search-time near food
csCmd = {'IPI.vrD', 'mm/IPI near Landmarks'; ... %R1        
         'abs(RS.vrV)', 'Speed near Landmarks'; ...
         'IPI.vrI*1000', 'IPI (ms) near Landmarks'};     
iZone = 2; %plot all active
plotAll({vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}, csCmd, iZone);

% search-time near food
csCmd = {'IPI.vrD', 'mm/IPI near Food'; ... %R1        
         'abs(RS.vrV)', 'Speed near Food'; ...
         'IPI.vrI*1000', 'IPI (ms) near Food'};     
iZone = 5; %plot all active
plotAll({vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}, csCmd, iZone);


%% build corr bar graphs

vsTrialPool = [vsTrialPool_E, vsTrialPool_L, vsTrialPool_P];

IPI_A = poolTrials_IPI(vsTrialPool, 1);
IPI_B = poolTrials_IPI(vsTrialPool, 2);
IPI_C = poolTrials_IPI(vsTrialPool, 3);
IPI_D = poolTrials_IPI(vsTrialPool, 4);

figure;

subplot 221; 
csZ = {IPI_A.vrCorrI_D, IPI_A.vrCorrV_D, IPI_A.vrCorrIV_D;...
       IPI_B.vrCorrI_D, IPI_B.vrCorrV_D, IPI_B.vrCorrIV_D;...
       IPI_C.vrCorrI_D, IPI_C.vrCorrV_D, IPI_C.vrCorrIV_D;...
       IPI_D.vrCorrI_D, IPI_D.vrCorrV_D, IPI_D.vrCorrIV_D};                 
plotBarMuSem(csZ, {'A', 'B', 'C', 'D'}, {'m', 'c', 'k'});
ylabel('Correlation');

subplot 222; 
csZ = {IPI_A.vrVafI_D, IPI_A.vrVafV_D, IPI_A.vrVafIV_D;...
       IPI_B.vrVafI_D, IPI_B.vrVafV_D, IPI_B.vrVafIV_D;...
       IPI_C.vrVafI_D, IPI_C.vrVafV_D, IPI_C.vrVafIV_D;...
       IPI_D.vrVafI_D, IPI_D.vrVafV_D, IPI_D.vrVafIV_D};     
plotBarMuSem(csZ, {'A', 'B', 'C', 'D'}, {'m', 'c', 'k'});
ylabel('VAF');

subplot 223; 
csZ = {IPI_A.vrCorrV_I, IPI_A.vrCorrA_I, IPI_A.vrCorrA_VHA;...
       IPI_B.vrCorrV_I, IPI_B.vrCorrA_I, IPI_B.vrCorrA_VHA;...
       IPI_C.vrCorrV_I, IPI_C.vrCorrA_I, IPI_C.vrCorrA_VHA;...
       IPI_D.vrCorrV_I, IPI_D.vrCorrA_I, IPI_D.vrCorrA_VHA};                 
plotBarMuSem(csZ, {'A', 'B', 'C', 'D'}, {'m', 'c', 'k'});
ylabel('Correlation');

subplot 224; 
csZ = {IPI_A.vrVafV_I, IPI_A.vrVafA_I, IPI_A.vrVafA_VHA;...
       IPI_B.vrVafV_I, IPI_B.vrVafA_I, IPI_B.vrVafA_VHA;...
       IPI_C.vrVafV_I, IPI_C.vrVafA_I, IPI_C.vrVafA_VHA;...
       IPI_D.vrVafV_I, IPI_D.vrVafA_I, IPI_D.vrVafA_VHA};         
plotBarMuSem(csZ, {'A', 'B', 'C', 'D'}, {'m', 'c', 'k'});
ylabel('VAF');



%% Dist per DIPI
% Dist per DIPI

csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csPlot = {'Early', 'Late', 'Probe'}; 
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csPlot = {'Stable', 'None', 'Unstable'}; 

csCmd = {'IPI.vrDsac', 'Esac/cm', '@(x)1./mean(x)', '@(x)0'; ...
         'IPI.vrDsac_LM', 'Esac/cm', '@(x)1./mean(x)', '@(x)0'; ...
         'IPI.vrDsac_F15', 'Esac/cm', '@(x)1./mean(x)', '@(x)0'};
plotAll(csTrials, csCmd, 1, csPlot);
suptitle('Active zone');
axis([.5 3.5 .05 .35]); set(gca, 'YTick', .05:.1:.35);
ylabel('');

%% Dist per DIPI
csCmd = {'IPI.vrDsac_LM', 'cm/Esac'; ...
         'IPI.vrIsac_LM', 'Esac interval (s)'};     
plotAll(csTrials, csCmd, [], csPlot);
suptitle('LM<3');

% Dist per DIPI
csCmd = {'IPI.vrDsac_F15', 'cm/Esac'; ...
         'IPI.vrIsac_F15', 'Esac interval (s)'};
plotAll(csTrials, csCmd, [], csPlot);
suptitle('F<15');

%% show traj
figure;
S = vsTrialPool_P(1);
imshow(imadjust(S.img0)); hold on;
mrX = S.mrX(1:5000,:);
mrY = S.mrY(1:5000,:);
for i=1:100:size(mrX,1)
    plot(mrX(i,[2,4]), mrY(i,[2,4]), 'b-');
    plot(mrX(i,[4,6]), mrY(i,[4,6]), 'r-');
    plot(mrX(i,[4]), mrY(i,[4]), 'k.');
end
% plot(mrX(:,2), mrY(:,2), 'r.'); hold on;
% plot(mrX(:,2), mrY(:,2), 'b-'); hold on;
plot(mrX(:,4), mrY(:,4), 'k:'); hold on; %midpoint

%% Sampling density

sem = @(x)std(x)/numel(x);
% csCmd = {'abs(IPI.vrD)', '<\Deltad>^{-1} (N/mm)'; ...
%          'abs(rad2deg(IPI.vrDHA))', '<\Delta\theta_H>^{-1} (N/deg)'; ...
%          'abs(rad2deg(IPI.vrDA))', '<\Delta\theta_T>^{-1} (N/deg)'; ...
%          'abs(RS.vrV)', 'Speed (cm/s)';...
%          'abs(rad2deg(RS.vrHAV))', '<\omega_H> (deg/s)';...
%          'abs(rad2deg(RS.vrV))', '<\theta_T> (deg/s)'};

csCmd = {'abs(IPI.vrD*10)', '<\Deltad> (mm/IPI)'; ...
         'abs(RS.vrV)', '<Speed> (cm/s)';...
         'abs(rad2deg(RS.vrAt))', '<\theta_T> (deg/s)'; ...
         'double(RS.vrV<0)', 'Prob. Back-swim'};

