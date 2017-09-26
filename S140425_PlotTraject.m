% S140425 plot trajectories
% start when entering the active zone

load D140330_Landmark;
load mlMask;

%%
angXaxis = -1.1590; %deg
rectCrop = [493 1083 312 902];

imgel = imadjust(imrotate(img0, angXaxis, 'nearest', 'crop'));
imgp = imadjust(imrotate(img0_P, angXaxis, 'nearest', 'crop'));


figure;

axis(rectCrop);

% draw trajectory
% limInt = stretchlim(img0(mlMask));
iPhase = 1;
nPhase = 3;
nTrial = 16;
iTrial = 1;
for iPhase=1:nPhase
    switch iPhase
        case 1, img = imgel;
        case 2, img = imgel;
        case 3, img = imgp;
    end
    for iTrial=1:nTrial
        subplot(nPhase, nTrial, iTrial + (iPhase-1)*nTrial);
        imshow(img); axis(rectCrop); hold on;
        csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P};

        vsTrials = csTrials{iPhase};
        S = poolTrials_RS(vsTrials(iTrial));
        vrX = S.vrX;
        vrY = S.vrY;

        %find the index when it gets to the food
        [~, ib] = min(S.vrDf);
        ia = find(~S.vlZone(1:ib), 1, 'last')-1;
        if isempty(ia), ia = 1; end
        viRng = ia:ib;
        plot(vrX(viRng), vrY(viRng), 'color', 'b');
        title(sprintf('P:%d, T:%d', iPhase, iTrial));
    end
end

%% pick three each and superimpose
miTrials = [3,6,9; 1,2,3; 2,5,15];
vcColor = 'rgb';
figure;
for iPhase=1:3
    vsTrials = csTrials{iPhase};
    switch iPhase
        case 1, img = imgel;
        case 2, img = imgel;
        case 3, img = imgp;
    end
    vi = miSelect(iPhase, :);
    subplot(1,3,iPhase);
    imshow(img); axis(rectCrop); hold on;
    for iPlot=1:3
        iTrial = miTrials(iPhase,iPlot);
        S = poolTrials_RS(vsTrials(iTrial));
        vrX = S.vrX;
        vrY = S.vrY;
        %find the index when it gets to the food
        [~, ib] = min(S.vrDf);
        ia = find(~S.vlZone(1:ib), 1, 'last')-1;
        if isempty(ia), ia = 1; end
        viRng = ia:ib;
        plot(vrX(viRng), vrY(viRng), 'color', vcColor(iPlot));
    end
end