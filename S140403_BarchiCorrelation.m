%% Correlation of trajectories
load D140330_Landmark;
load D140330_Landmark_Probe;

%% compute visit count and correlation
mlMask = getImageMask(img0, [0 50], 'CENTRE');

viAnimal_E = poolVecFromStruct(vsTrialPool_E, 'iAnimal');
viAnimal_L = poolVecFromStruct(vsTrialPool_L, 'iAnimal');
cmVD_Ea = {}; %visit density
cmVD_Eb = {};
cmVD_La = {};
cmVD_Lb = {};
for iAnimal=1:4
    viTrial_E = find(viAnimal_E == iAnimal);
    viTrial_L = find(viAnimal_L == iAnimal);    
    viTrial_Ea = viTrial_E(1:round(end/2));
    viTrial_Eb = viTrial_E(round(end/2)+1:end);
    viTrial_La = viTrial_L(1:round(end/2));
    viTrial_Lb = viTrial_L(round(end/2)+1:end);
    
    [~, cmVD_Ea{iAnimal}] = calcVisitCount(vsTrialPool_E(viTrial_Ea), img0, mlMask);
    [~, cmVD_Eb{iAnimal}] = calcVisitCount(vsTrialPool_E(viTrial_Eb), img0, mlMask);
    [~, cmVD_La{iAnimal}] = calcVisitCount(vsTrialPool_L(viTrial_La), img0, mlMask);
    [~, cmVD_Lb{iAnimal}] = calcVisitCount(vsTrialPool_L(viTrial_Lb), img0, mlMask);
end

%% find max cross-correlation score per animal
threshXA = .6;

vrXC_E = zeros(4,1); %max xcorr2
vrXC_L = zeros(4,1);
vrXA_E = zeros(4,1); %area above 60%
vrXA_L = zeros(4,1);
mrXC = zeros(4,4);  % points above threshold over days
mrXA = zeros(4,4);

for iAnimal=1:4
    [vrXC_E(iAnimal), vrXA_E(iAnimal)] = ...
        calcXcorr2Max(cmVD_Ea{iAnimal}, cmVD_Eb{iAnimal}, threshXA);
    
    [vrXC_L(iAnimal), vrXA_L(iAnimal)] = ...
        calcXcorr2Max(cmVD_La{iAnimal}, cmVD_Lb{iAnimal}, threshXA);

    [mrXC(1, iAnimal), mrXA(1, iAnimal)] = calcXcorr2Max(cmVD_Ea{iAnimal}, cmVD_Lb{iAnimal}, threshXA);
    [mrXC(2, iAnimal), mrXA(2, iAnimal)] = calcXcorr2Max(cmVD_Eb{iAnimal}, cmVD_Lb{iAnimal}, threshXA);
    [mrXC(3, iAnimal), mrXA(3, iAnimal)] = calcXcorr2Max(cmVD_La{iAnimal}, cmVD_Lb{iAnimal}, threshXA);
    [mrXC(4, iAnimal), mrXA(4, iAnimal)] = calcXcorr2Max(cmVD_Lb{iAnimal}, cmVD_Lb{iAnimal}, threshXA);
end
figure;
subplot 221;
bar(1:4, [vrXC_E, vrXC_L], 1);
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D'});
xlabel('Animal');
ylabel('Maximum correlation');
title('(Ea x Eb) vs. (La x Lb)');
set(gca, 'XLim', [.5 4.5]);

subplot 222;
bar(1:4, [vrXA_E, vrXA_L], 1);
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D'});
xlabel('Animal');
ylabel('Points above threshold');
title('(Ea x Eb) vs. (La x Lb)');
set(gca, 'XLim', [.5 4.5]);

subplot 223;
plot(mrXC, ':.');
set(gca, 'XTickLabel', {'Ea', 'Eb', 'La', 'Lb'});
xlabel('Time');
ylabel('Maximum correlation');
title('Correlation with Lb');
set(gca, 'XLim', [.5 4.5]);

subplot 224;
plot(mrXA, ':.');
set(gca, 'XTickLabel', {'Ea', 'Eb', 'La', 'Lb'});
xlabel('Time');
ylabel('Points above threshold');
title('Correlation with Lb');
set(gca, 'XLim', [.5 4.5]);

%% show the trajectories
mlMask1 = imresize(mlMask, 1/25);
figure; 
for iAnimal=1:4
    offset = (iAnimal-1)*4;
    subplot(4,4,1 + offset);
    imagesc(imCropMask(cmVD_Ea{iAnimal}, mlMask1));
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
    
    subplot(4,4,2 + offset);
    imagesc(imCropMask(cmVD_Eb{iAnimal}, mlMask1));
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
    
    subplot(4,4,3 + offset);
    imagesc(imCropMask(cmVD_La{iAnimal}, mlMask1));
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
    
    subplot(4,4,4 + offset);
    imagesc(imCropMask(cmVD_Lb{iAnimal}, mlMask1)); 
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
end
suptitle('Visit Density');

%% show Autocorrelation
mlMask1 = imresize(mlMask, 1/25);
figure; 
for iAnimal=1:4
    offset = (iAnimal-1)*4;
    subplot(4,4,1 + offset);
    imagesc(xcorr2(imCropMask(cmVD_Ea{iAnimal}, mlMask1)));
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
    
    subplot(4,4,2 + offset);
    imagesc(xcorr2(imCropMask(cmVD_Eb{iAnimal}, mlMask1)));
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
    
    subplot(4,4,3 + offset);
    imagesc(xcorr2(imCropMask(cmVD_La{iAnimal}, mlMask1)));
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
    
    subplot(4,4,4 + offset);
    imagesc(xcorr2(imCropMask(cmVD_Lb{iAnimal}, mlMask1))); 
    axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
end
suptitle('XCorr2 Visit Density');