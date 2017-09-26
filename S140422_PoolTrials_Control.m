%%-----------------------------------------
% IMPORT data and pool
% Pooled-trial analysis


for iPhase=1:3
    switch iPhase
        case 1
            dname = 'C:\expr\Control_Empty\'; %run probe first and go onto 
        case 2
            dname = 'C:\expr\Control_Rand\';
        case 3
            dname = 'C:\expr\Control_RandWide\';
    end
    vsTrials = [];
    warning off;

    % Collect trials
    vsFnames = dir([dname '*_Track.mat']);
    nFiles = numel(vsFnames);
    for i=1:nFiles
        fname = vsFnames(i).name;
        disp(fname);
        vsTrials = [vsTrials, importTrial(load([dname fname]))];
    end

    switch iPhase
        case 1
            vsTrialPool_NP = vsTrials;
        case 2
            vsTrialPool_RP = vsTrials;
        case 3
            vsTrialPool_WP = vsTrials;
    end
end

save D140422_Controls vsTrialPool_NP vsTrialPool_RP vsTrialPool_WP;

%% check for image alignments
xy0 = [];
rectCrop = [493 1083 312 902];
rect1 = [rectCrop(1), rectCrop(3), diff(rectCrop(1:2)), diff(rectCrop(3:4))];
for iPhase=1:3
    disp(iPhase);
    switch iPhase
        case 1
            vsTrials = vsTrialPool_NP;            
        case 2
            vsTrials = vsTrialPool_RP;
        case 3
            vsTrials = vsTrialPool_WP;
    end
    xy0 = vsTrials(1).xy0;
    
    for iTrial=1:numel(vsTrials)
        S = vsTrials(iTrial);
        hold off;
        imshow(imadjust(S.img0));        
        hold on; plot(xy0(1), xy0(2), 'r.');
        title(sprintf('Phase: %d, Trial: %d, %s', iPhase, iTrial, S.dataID));
        imrect(gca, rect1);
        pause;
    end
end

