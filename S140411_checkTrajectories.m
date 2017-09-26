% S140411_checkTrajectories

%% load individual trials and check trajectory

dname = 'C:\expr\Landmark_Probe\'; %run probe first and go onto 

% dname = 'C:\expr\Landmark_Learning\';

% Collect trials
vsFnames = dir([dname '*_Track.mat']);
nFiles = numel(vsFnames) %16 files

iFile=0;
%% Kinematics calculation
figure;
warning off;


for iFile=1:nFiles
% iFile = iFile+1; %out of 16
clf;
fname = vsFnames(iFile).name;
S = importTrial(load([dname fname]));

subplot 121; 
imshow(imadjust(S.img0)); hold on;
plot(S.XH, S.YH);
plot(S.xy0(1), S.xy0(2), 'r.');
plot(S.XH(1), S.YH(1), 'r*');
plot(S.XH(end), S.YH(end), 'g*');
title(S.dataID);

subplot 122; 
imshow(imadjust(S.img0)); hold on;
plot(S.vrX, S.vrY);
plot(S.xy0(1), S.xy0(2), 'r.');
plot(S.vrX(1), S.vrY(1), 'r*');
plot(S.vrX(end), S.vrY(end), 'g*');
title('rotated');
pause(.5);
end

%%

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