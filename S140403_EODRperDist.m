%% Correlation of trajectories
load D140330_Landmark;

% Variables in the landmark struct
%   TEOD, EODR, EODA, EODRq, EODRz, vtEOD, ...
%   XH, YH, VEL, ACC, ANG, AVEL, ...
%   viESAC, vrESAC, vtESAC, vxESAC, vyESAC, img0, ...
%   dataID, pathLen_cm, xyFood, xy0, xyStart, ...
%   vlZone, Rfood, Rcentre, vrX, vrY, iAnimal, duration);

%
get_E = @(x)poolVecFromStruct(vsTrialPool_E, x);
get_L = @(x)poolVecFromStruct(vsTrialPool_L, x);
get_P = @(x)poolVecFromStruct(vsTrialPool_P, x);
region_E = @(ml)getRegion(vsTrialPool_E, ml);
region_L = @(ml)getRegion(vsTrialPool_L, ml);
region_P = @(ml)getRegion(vsTrialPool_P, ml);

bootCV = @(x)[std(x)/mean(x); bootci(1000, {@(y)std(y)/mean(y), x})];
bootSD = @(x)[std(x); bootci(1000, {@(y)std(y), x})];
bootMean = @(x)[mean(x); bootci(1000, {@(y)mean(y), x}, 'type', 'cper')];

%% compute visit count and correlation
strRegion = 'CENTRE';
limDist = [0 50];

mlMask = getImageMask(img0, limDist, strRegion);

EODPerDist = @(x)poolVecFromStruct(x, 'EODR') ./ abs(poolVecFromStruct(x, 'VEL'));

vrNpd_E = EODPerDist(vsTrialPool_E);
vrNpd_L = EODPerDist(vsTrialPool_L);
vrNpd_P = EODPerDist(vsTrialPool_P);

vrNpd_E1 = vrNpd_E(region_E(mlMask));
vrNpd_L1 = vrNpd_L(region_L(mlMask));
vrNpd_P1 = vrNpd_P(region_P(mlMask));

figure; hold on;
xi = 0:.01:6;
ksdensity(vrNpd_E1, xi, 'function', 'survivor'); 
ksdensity(vrNpd_L1, xi, 'function', 'survivor');
ksdensity(vrNpd_P1, xi, 'function', 'survivor');
h = get(gca, 'Children'); 
set(h(3), 'Color', 'r');
set(h(2), 'Color', 'b');
set(h(1), 'Color', 'g');
legend({'E', 'L', 'P'});
xlabel('EOD #/Dist (#/cm)');
title(sprintf('%s, %d~%d cm', strRegion, limDist(1), limDist(2)));


%% plot swim speed per animal
% bootMean = @(x)[mean(x); bootci(1000, {@(y)mean(y), x}, 'type', 'cper')];

pixpercm = 1053.28/(sqrt(2)*100);
strRegion = 'CENTRE';
limDist = [0 50];
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

mcSpeed = cell(4, 3);
for iAnimal=1:4
    mcSpeed{iAnimal, 1} = abs(getVecFromAnimal(vsTrialPool_E, 'VEL', iAnimal, mlMask));
    mcSpeed{iAnimal, 2} = abs(getVecFromAnimal(vsTrialPool_L, 'VEL', iAnimal, mlMask));
    mcSpeed{iAnimal, 3} = abs(getVecFromAnimal(vsTrialPool_P, 'VEL', iAnimal, mlMask));   
end

figure; 
plotCellErrorbar(mcSpeed, 'VEL', @(x)x/pixpercm);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
ylabel('Swim speed (cm/s)');
title(regionStr);
set(gca, 'YLim', [0 15])

set(gca, 'XTickLabel', {'A', 'B', 'C', 'D'});
title(strRegion);
ylabel('Swim speed (cm/s)');


%% ---------------------------------------------------------
% show swim speed per phase
pixpercm = 1053.28/(sqrt(2)*100);

strRegion = 'CENTRE';
limDist = [0 50];
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

csSwimSpeedAll = cell(3,1);
csSwimSpeedAll{1} = abs(poolVecFromStruct(vsTrialPool_E, 'VEL', mlMask));
csSwimSpeedAll{2} = abs(poolVecFromStruct(vsTrialPool_L, 'VEL', mlMask));
csSwimSpeedAll{3} = abs(poolVecFromStruct(vsTrialPool_P, 'VEL', mlMask));

figure; 
plotCellErrorbar(csSwimSpeedAll, 'VEL', @(x)x/pixpercm);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
ylabel('Swim speed (cm/s)');
title(regionStr);
set(gca, 'YLim', [0 15])


%% divide survivor function to forward abnd backward

%%---------------------------------------------------------
%%
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


%%---------------------------------------------------
%% # samples / cm

Strial = vsTrialPool_E(1); % Trial to analyze
strPhase = 'Early A1 ';

midVal = @(x)(x(1:end-1) + x(2:end))/2;

vtEOD = poolVecFromStruct(Strial, 'vtEOD');
vrX = poolVecFromStruct(Strial, 'XH');
vrY = poolVecFromStruct(Strial, 'YH');
vrT = poolVecFromStruct(Strial, 'TEOD');
vrV = poolVecFromStruct(Strial, 'VEL');

% interpolate the locatoin of the EOD pulse
vrXe = interp1(vrT, vrX, vtEOD, 'spline');
vrYe = interp1(vrT, vrY, vtEOD, 'spline');
vrD = sqrt(diff(vrXe).^2 + diff(vrYe).^2);
vrTd = midVal(vtEOD);
vrXd = midVal(vrXe);
vrYd = midVal(vrYe);
vrId = diff(vtEOD);
vrVd = interp1(vrT, vrV, vrTd, 'spline');

%% visualize
figure; plot(vrT, vrX, '-');
% interpolate the locatoin of pulse
hold on; plot(vtEOD, vrXe, 'r.');


figure; imshow(imadjust(img0)); hold on;
plot(vrXd(1:10:end), vrYd(1:10:end), '.');

%corr score: -0.2780 
figure; plot(vrVd, 1./vrId, '.'); xlabel('vel'); ylabel('inst R'); 
figure; plot(abs(vrVd), zscore(1./vrId), '.'); xlabel('Speed (cm/s)'); ylabel('Inst. EODRz'); 
figure; plot(abs(vrVd), vrD, '.'); xlabel('Speed (cm/s)'); ylabel('Dist per IPI (cm)'); 
figure; plot(vrD, zscore(1./vrId), '.'); xlabel('Speed (cm/s)'); ylabel('Inst. EODRz'); 

%% location-dependent assay of the dist per trajectory
% 

strRegion = 'CENTRE';
limDist = [0 50];
regionStr = sprintf('within %d~%d of %s', limDist(1), limDist(2), strRegion);
mlMask = getImageMask(img0, limDist, strRegion);

vlRegion = mlMask(sub2ind(size(mlMask), round(vrYd), round(vrXd)));

% filter region do this at last
% vtEOD = vtEOD(viRegion);
vrXd1 = vrXd(vlRegion);
vrYd1 = vrYd(vlRegion);
vrId1 = vrId(vlRegion);
vrD1 = vrD(vlRegion);
vrVd1 = vrVd(vlRegion);

figure; ksdensity(vrD1);
title([strPhase, ', ', regionStr]);
xlabel('Dist per IPI (cm)');

%%
% strVar = 'vrI';   f1=@(x)diff(x)*1000
strVar = 'vrD'; f1=@(x)x;
strRegion = 'CENTRE';
limDist = [0 50];
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

csDistIPI = cell(4,3);
for iAnimal=1:4
    csDistIPI{iAnimal,1} = poolDistPerIPI(vsTrialPool_E, iAnimal, mlMask);
    csDistIPI{iAnimal,2} = poolDistPerIPI(vsTrialPool_L, iAnimal, mlMask);
    csDistIPI{iAnimal,3} = poolDistPerIPI(vsTrialPool_P, iAnimal, mlMask);
end

figure; 
subplot 121; %show each animals
plotCellErrorbar(csDistIPI, strVar, @(x)f1(x));
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D'});
ylabel('IPI per dist. (1/cm)');   set(gca, 'YLim', [0 .25]);
% ylabel('IPI (ms)'); set(gca, {'YLim', 'YTick'}, {[13 15], 13:.5:15});
title(regionStr);


subplot 122;
csDistIPIall = cell(3,1);
csDistIPIall{1} = poolDistPerIPI(vsTrialPool_E, [], mlMask);
csDistIPIall{2} = poolDistPerIPI(vsTrialPool_L, [], mlMask);
csDistIPIall{3} = poolDistPerIPI(vsTrialPool_P, [], mlMask);
plotCellErrorbar(csDistIPIall, strVar, @(x)f1(x));
set(gca, 'XTickLabel', {'E', 'L', 'P'});
ylabel('IPI per dist. (1/cm)');  set(gca, 'YLim', [0 .25]);
% ylabel('IPI (ms)'); set(gca, {'YLim', 'YTick'}, {[13 15], 13:.5:15});
title(regionStr);

% figure; plotXcorr(csDistIPI{iAnimal,1}.vrD);
% calcCorrTau(csDistIPI{iAnimal,1}.vrD)

%% significance test
[h p] = kstest2_jjj(csDistIPIall{1}.vrD, csDistIPIall{2}.vrD)
[h p] = kstest2_jjj(csDistIPIall{2}.vrD, csDistIPIall{3}.vrD)

a=poolVecFromStruct(csDistIPI_E{2}, 'vrDI');
calcDistAsym(a)
a=poolVecFromStruct(csDistIPI_L{2}, 'vrDI');
calcDistAsym(a)
a=poolVecFromStruct(csDistIPI_P{2}, 'vrDI');
calcDistAsym(a)

a=poolVecFromStruct(csDistIPI_E, 'vrDI');

strRegion = 'CENTRE';
limDist = [0 50];
mlMask = getImageMask(img0, limDist, strRegion);

vrVec = getVecFromAnimal(csDistIPI_E, 'vrDI', 1:4, mlMask);
calcDistAsym(vrVec)


%% negative vs positive velocity and distance

strRegion = 'CENTRE';
limDist = [0 50];
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

csDistIPIall_E = poolDistPerIPI(vsTrialPool_E, [], mlMask);
vrD_E = poolVecFromStruct(csDistIPIall_E, 'vrD');
vrV_E = poolVecFromStruct(csDistIPIall_E, 'vrV');

csDistIPIall_L = poolDistPerIPI(vsTrialPool_L, [], mlMask);
vrD_L = poolVecFromStruct(csDistIPIall_L, 'vrD');
vrV_L = poolVecFromStruct(csDistIPIall_L, 'vrV');

csDistIPIall_P = poolDistPerIPI(vsTrialPool_P, [], mlMask);
vrD_P = poolVecFromStruct(csDistIPIall_P, 'vrD');
vrV_P = poolVecFromStruct(csDistIPIall_P, 'vrV');

cvV = cell(3,2);
cvV{1,1} = vrD_E(vrV_E>0);
cvV{1,2} = vrD_E(vrV_E<0);
cvV{2,1} = vrD_L(vrV_L>0);
cvV{2,2} = vrD_L(vrV_L<0);
cvV{3,1} = vrD_P(vrV_P>0);
cvV{3,2} = vrD_P(vrV_P<0);

figure; plotCellErrorbar(cvV);
title(['Forward (blue) vs. Backward (red) swimming, ' regionStr]);
set(gca, 'XTickLabel', {'E', 'L', 'P'});
ylabel('Dist per IPI (cm)');
set(gca, 'YLim', [0 2]);