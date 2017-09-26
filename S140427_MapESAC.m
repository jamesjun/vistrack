% S140527_Map ESAC

rectCrop = [493 1083 312 902];

titleStr = 'IPI/mm, Late trials';

csPhase = {'Early', 'Late', 'Probe'};
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P};

img0 = [];
figure; 
for iPhase = 1:3
    subplot(1,3,iPhase);
    vsTrial = csTrials{iPhase};
%     S = poolDIPIperDist(vsTrial, []); %replaced by poolTrialsIPI
    S = poolTrials_IPI(vsTrial);
    [mnVisit, mnVisit1, mnI] = mapEscan(S, img0);
    imshow(mnI);
    axis(rectCrop);
    title(csPhase(iPhase));
end
suptitle('<mm/IPI>^-1');