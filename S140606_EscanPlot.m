%% analyze down and up-state shot noise
S = load('D:\Expr\fishC\Cspont_120608_Ts');
vrT = S.V120608_freeswim_000_Ch1.times;
vrI = differentiate3(vrT);
vrDI = differentiate3(vrI);

[vi, thresh, vrAmpl, viAmpl] = findDIsac(vrDI);

figure; 
plot(vrT, vrDI * 1000, 'b', vrT(vi), vrDI(vi) * 1000, 'r.');
ylabel('\DeltaIPI (ms)');
xlabel('Time (s)');


vrEscanInterval = diff(vrT(vi));
figure; ksdensity(vrEscanInterval, 'function', 'pdf'); 
xlabel('Time (s)'); 
title('E-scan interval');
set(gca, 'XLim', [0 5]);
% set(gca, 'YScale', 'log');
set(gca, 'XLim', [0 5]);

figure; ksdensity(vrEscanInterval, 'function', 'pdf'); 

%% Plot distribution
thresh = -40e-6; %in sec, common threshold +/- 5 us
TLIM_dn = [2.74, 2.78] * 10000;
TLIM_up = [2.785, 2.825] * 10000;

vlUp = vrT >= TLIM_up(1) & vrT < TLIM_up(2);
vlDn = vrT >= TLIM_dn(1) & vrT < TLIM_dn(2);
vlAll = vrT >= TLIM_dn(1) & vrT < TLIM_up(2);

[vi_up, thresh_up, vrAmpl_up, vnDur_up] = findDIsac(vrDI(vlUp), thresh*2);
[vi_dn, thresh_dn, vrAmpl_dn, vnDur_dn] = findDIsac(vrDI(vlDn), thresh*2);
disp(thresh_up)
disp(thresh_dn)

% plot up and down states
figure; AX = [];
subplot 211; AX(1) = gca;
plot(vrT(vlAll), vrI(vlAll)*1000, 'k', vrT(vlUp), vrI(vlUp)*1000, 'r', vrT(vlDn), vrI(vlDn)*1000, 'b');
xlabel('Time');
ylabel('IPI (ms)');
axis([2.74e4 2.825e4 13 21])
subplot 212; AX(2) = gca;
plot(vrT(vlAll), vrDI(vlAll)*1000, 'k', vrT(vlUp), vrDI(vlUp)*1000, 'r', vrT(vlDn), vrDI(vlDn)*1000, 'b');
xlabel('Time');
ylabel('\DeltaIPI (ms)');
linkaxes(AX, 'x');
% axis([2.74e4 2.825e4 13 21])

xi = 0:.001:10;
figure; 
vrIes_up = diff(vrT(vi_up));
ksdensity(vrIes_up, xi, 'function', 'survivor'); 
xlabel('E-scan Interval (s)'); ylabel('Survivor');
axis([0 2 1e-3 10]); set(gca, 'YScale', 'log');
title('Up-state survivor function');

figure; 
vrIes_dn = diff(vrT(vi_dn));
ksdensity(vrIes_dn, xi, 'function', 'survivor'); 
xlabel('Escan Interval Dn state'); ylabel('Survivor');
axis([0 2 1e-3 10]); set(gca, 'YScale', 'log');
title('Down-state survivor function');

%% Fit distribution
trefrac = .41;
tdur = 1.5;

UP = vrIes_up - trefrac;
UP(UP<0 | UP > tdur) = [];
dfittool(UP);
%mu        0.314461  0.0172066

DN = vrIes_dn - trefrac;
DN(DN<0 | DN > tdur) = [];
dfittool(DN);
% mu       0.233786  0.0146403

%% Manual fit
% [f]  = ksdensity(diff(vrT(vi_up)), xi, 'function', 'survivor'); 
% xi = 0:.001:3;
% vrIes_up = diff(vrT(vi_up));
% vrIes_dn = diff(vrT(vi_dn));
% 
% f_up = ksdensity(vrIes_up, xi);
% f_dn = ksdensity(vrIes_dn, xi);
% s_up = 1-cumsum(f_up)/sum(f_up);
% s_dn = 1-cumsum(f_dn)/sum(f_dn);
% 
% figure; plot(xi, log(f_up), 'r', xi, log(f_dn), 'b');

%% Fig. 8D plot. E-scan density before vs. after
cStrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csXlabel = {'Early', 'Late', 'Probe'};        
iZone = 1;
iAnimal = [];
cmDsac_Az = {};
cmDsac_LM = {};
cmDsac_F15 = {};

for iPhase = 1:numel(cStrials)
    vsTrials = cStrials{iPhase};
    for iTrial = 1:numel(vsTrials)
        IPI = poolTrials_IPI(vsTrials(iTrial));
        cmDsac_Az{iTrial, iPhase} = 1./nanmean(IPI.vrDsac_Az);
        cmDsac_LM{iTrial, iPhase} = 1./nanmean(IPI.vrDsac_LM);
        cmDsac_F15{iTrial, iPhase} = 1./nanmean(IPI.vrDsac_F15);
    end
end


trDsac_Az = mean_bootsem_cm(cmDsac_Az);
trDsac_LM = mean_bootsem_cm(cmDsac_LM);
trDsac_F15 = mean_bootsem_cm(cmDsac_F15);

figure; mean_bootsem_cm(cmDsac_LM); title('Near Landmark'); ylim([0 .6]);
figure; mean_bootsem_cm(cmDsac_F15); title('Near food'); ylim([0 .6]);

%%
figure;
subplot 131;
plotErrorbar(trDsac_Az); 
set(gca, 'XTickLabel', csXlabel);
ylabel('E-sac density (cm^{-1})');
title('Active zone')

    
subplot 132;
plotErrorbar(trDsac_F15); 
set(gca, 'XTickLabel', csXlabel);
ylabel('E-sac density (cm^{-1})');
title('Around food');

subplot 133;
plotErrorbar(trDsac_LM); 
set(gca, 'XTickLabel', csXlabel);
ylabel('E-sac density (cm^{-1})');
title('LM');