% Fig. 7 and S3 combined. backward swim plot

%% load
load D140330_Landmark;
cStrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csTrials = {'Early', 'Late', 'Probe'}; 

%% process
csXlabel = {'Early', 'Late', 'Probe'};        
iZone = 1;
iAnimal = [];
clear trFB*;
for iPhase = 1:numel(cStrials)
    vsTrials = cStrials{iPhase};
    clear mrFB*;
    for iTrial = 1:numel(vsTrials)
        IPI = poolTrials_IPI(vsTrials(iTrial), iAnimal, [], iZone);
        vrV = abs(IPI.vrV);
        vrR = 1./IPI.vrI;
        vrSd = 1./IPI.vrD;    
        mrFB_Sd(iTrial,1) = 1./nanmean(IPI.vrD(IPI.vrV>0 & IPI.vlZ0));
        mrFB_Sd(iTrial,2) = 1./nanmean(IPI.vrD(IPI.vrV<0 & IPI.vlZ0));
        mrFB_R(iTrial,1) = nanmean(vrR(IPI.vrV>0 & IPI.vlZ0));
        mrFB_R(iTrial,2) = nanmean(vrR(IPI.vrV<0 & IPI.vlZ0));
        mrFB_V(iTrial,1) = nanmean(vrV(IPI.vrV>0 & IPI.vlZ0));
        mrFB_V(iTrial,2) = nanmean(vrV(IPI.vrV<0 & IPI.vlZ0));
        %mrFB_corr_VSd(iTrial,1) = calcVAF(vrV, vrSd, IPI.vlZ0 & IPI.vrV>0);         
        %mrFB_corr_VSd(iTrial,2) = calcVAF(vrV, vrSd, IPI.vlZ0 & IPI.vrV<0); 
        %mrFB_corr_RSd(iTrial,1) = calcVAF(vrR, vrSd, IPI.vlZ0 & IPI.vrV>0);
        %mrFB_corr_RSd(iTrial,2) = calcVAF(vrR, vrSd, IPI.vlZ0 & IPI.vrV<0);
        %mrFB_corr_RV(iTrial,1) = calcVAF(vrR, vrV, IPI.vlZ0 & IPI.vrV>0);
        %mrFB_corr_RV(iTrial,2) = calcVAF(vrR, vrV, IPI.vlZ0 & IPI.vrV<0);
    end    
    trFB_Sd(:,:,iPhase) = mean_bootsem(mrFB_Sd);
    trFB_R(:,:,iPhase) = mean_bootsem(mrFB_R);
    trFB_V(:,:,iPhase) = mean_bootsem(mrFB_V);
    %trFB_corr_VSd(:,:,iPhase) = mean_bootsem(mrFB_corr_VSd);
    %trFB_corr_RSd(:,:,iPhase) = mean_bootsem(mrFB_corr_RSd);
    %trFB_corr_RV(:,:,iPhase) = mean_bootsem(mrFB_corr_RV);
    
    kstest2_disp(mrFB_Sd(:,1), mrFB_Sd(:,2), sprintf('Phase%d, SD', iPhase));
    kstest2_disp(mrFB_R(:,1), mrFB_R(:,2), sprintf('Phase%d, R', iPhase));
    kstest2_disp(mrFB_V(:,1), mrFB_V(:,2), sprintf('Phase%d, V', iPhase));
end

trFB_Sd = permute(trFB_Sd, [3,1,2]);
trFB_R = permute(trFB_R, [3,1,2]);
trFB_V = permute(trFB_V, [3,1,2]);
%trFB_corr_VSd = permute(trFB_corr_VSd, [3,1,2]);
%trFB_corr_RSd = permute(trFB_corr_RSd, [3,1,2]);
%trFB_corr_RV = permute(trFB_corr_RV, [3,1,2]);

% Plot
csXlabel = {'Early', 'Late', 'Probe'};

%%
figure;
subplot 131;
plotErrorbar(trFB_Sd); 
set(gca, 'XTickLabel', csXlabel);
ylabel('Sampling density (EOD/cm)');
set(gca, 'YTick', 0:10:40);
ylim([0,40]);

subplot 132;
plotErrorbar(trFB_R); 
set(gca, 'XTickLabel', csXlabel);
ylabel('EOD Rate (Hz)');
ylim([70 78]); 
set(gca, 'YTick', 70:2:78);

subplot 133;
plotErrorbar(trFB_V); 
set(gca, 'XTickLabel', csXlabel);
ylabel('Speed (cm/s)');
set(gca, 'YTick', 0:4:16);
ylim([0,16]);

set(gcf, 'Color', 'w');
colormap hot;


%%
subplot 234;
plotErrorbar(trFB_corr_VSd); 
set(gca, 'XTickLabel', csXlabel);
title('VAF V-Sd');

subplot 235;
plotErrorbar(trFB_corr_RSd); 
set(gca, 'XTickLabel', csXlabel);
title('VAF R-Sd');

subplot 236;
plotErrorbar(trFB_corr_RV); 
set(gca, 'XTickLabel', csXlabel);
title('VAF R-V');

set(gcf, 'Color', 'w');