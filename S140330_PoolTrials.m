%%-----------------------------------------
% IMPORT data and pool
% Pooled-trial analysis


for iPhase=1:2
    if iPhase==1
        dname = 'C:\expr\Landmark_Probe\'; %run probe first and go onto 
    else
        dname = 'C:\expr\Landmark_Learning\';
    end
    vsTrials = [];
    warning off;

    % Collect trials
    vsFnames = dir([dname '*_Track.mat']);
    nFiles = numel(vsFnames);
    for i=1:nFiles
        fname = vsFnames(i).name;
        disp(fname);
        S = importTrial(load([dname fname]));
        vsTrials = [vsTrials, S];
    end

    %get the first image
    viAnimal = poolVecFromStruct(vsTrials, 'iAnimal');
    iTrialA1 = find(viAnimal==1, 1, 'first');
    img0 = vsTrials(iTrialA1).img0;
    img0a = imadjust(img0);

    figure; imshow(img0a); title('Map for A');

    figure; suptitle('Learning curve');
    for iAnimal=1:4
        subplot(2,2,iAnimal);
        vrDist = poolVecFromStruct(vsTrials(viAnimal==iAnimal), 'pathLen_cm');
        bar(log10(vrDist));
        title(sprintf('Animal %c', 'A'+iAnimal-1));    
        ylabel('log10 Distance (cm)');
        xlabel('Trial #');    
    end

    if iPhase==1
        vsTrialPool_P = vsTrials; %select all probe trials
        img0_P = img0; 
        vrDist_P = poolVecFromStruct(vsTrialPool_P, 'pathLen_cm');
        vrDur_P = poolVecFromStruct(vsTrialPool_P, 'duration');
    end
end
% save D140330_Landmark_Probe vsTrialPool_P img0;

% Seperate early vs. late learning trials
viSession_E = [1:2];
viSession_L = [7:12];
quantLim = [1/8 7/8];

viTrials_E = (viSession_E([1 end]) + [-1 0]) * 4 + [1 0];
viTrials_L = (viSession_L([1 end]) + [-1 0]) * 4 + [1 0];
viTrials_E = viTrials_E(1):viTrials_E(end);
viTrials_L = viTrials_L(1):viTrials_L(end);

vsTrialAll_E = [];
vsTrialAll_L = [];
viAnimal = poolVecFromStruct(vsTrials, 'iAnimal');
for iAnimal=1:4
    S = vsTrials(viAnimal == iAnimal);
    vsTrialAll_E = [vsTrialAll_E, S(viTrials_E)];
    vsTrialAll_L = [vsTrialAll_L, S(viTrials_L)];
end

vrDist_E = poolVecFromStruct(vsTrialAll_E, 'pathLen_cm');
limDist_E = quantile(vrDist_E, quantLim);
vlTrial_E = vrDist_E >= limDist_E(1) & vrDist_E <= limDist_E(2);

vrDist_L = poolVecFromStruct(vsTrialAll_L, 'pathLen_cm');
limDist_L = quantile(vrDist_L, quantLim);
vlTrial_L = vrDist_L >= limDist_L(1) & vrDist_L <= limDist_L(2);

vrDur_E = poolVecFromStruct(vsTrialAll_E, 'duration');
vrDur_L = poolVecFromStruct(vsTrialAll_L, 'duration');

vsTrialPool_E = vsTrialAll_E(vlTrial_E);
vsTrialPool_L = vsTrialAll_L(vlTrial_L);

save D140330_Landmark vsTrialPool_E vsTrialPool_L vsTrialPool_P vrDur_E vrDur_L vrDur_P vrDist_E vrDist_L vrDist_P img0 img0_P;

%% fraction of trials
bootCV = @(x)[std(x)/mean(x); bootci(1000, {@(y)std(y)/mean(y), x})];
bootSD = @(x)[std(x); bootci(1000, {@(y)std(y), x})];
bootMean = @(x)[mean(x); bootci(1000, {@(y)mean(y), x})];
dispMuSd = @(vec)disp([vec(1), diff(vec([2 3]))/2]);

% Fractino of trials over the limit
dispMuSd(bootMean(vrDist_E > 1000)) %0.4063 +/- 0.1719
dispMuSd(bootMean(vrDist_L > 1000)) %0.0521 +/- 0.0469
dispMuSd(bootMean(vrDur_E > 100)) %0.3750 +/- 0.1719
dispMuSd(bootMean(vrDur_L > 100)) %0.0417 +/- 0.0469

%note: this doesn't capture aborted trials. timed-out trials needs to be
%stated

%% Fig. 2. compare e-saccades distribution between early vs. late learning

vrEODRz_E = poolVecFromStruct(vsTrialPool_E, 'EODRz');
vrEODRz_L = poolVecFromStruct(vsTrialPool_L, 'EODRz');
vrEODA_E = poolVecFromStruct(vsTrialPool_E, 'EODA');
vrEODA_L = poolVecFromStruct(vsTrialPool_L, 'EODA');

figure; suptitle('Survivor log10 EODA');
hold on;
ksdensity(log10(vrEODA_E(vrEODA_E>0)), 'function', 'survivor');
ksdensity(log10(abs(vrEODA_E(vrEODA_E<0))), 'function', 'survivor');
ksdensity(log10(vrEODA_L(vrEODA_L>0)), 'function', 'survivor');
ksdensity(log10(abs(vrEODA_L(vrEODA_L<0))), 'function', 'survivor');
a = get(gca, 'Children'); 
set(a(1), 'Color', 'r');
set(a(2), 'Color', 'm');
set(a(3), 'Color', 'b');
set(a(4), 'Color', 'c');
legend({'EODA+ E', 'EODA- E', 'EODA+ L', 'EODA- L'});

figure; suptitle('Survivor EODA Late');
qqplot(log10(vrEODA_L(vrEODA_L>0)), log10(abs(vrEODA_L(vrEODA_L<0)))); 
xlabel('log10 EODA+ L'); ylabel('log10 EODA- L');

