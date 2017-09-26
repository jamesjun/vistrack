%% Grid based amp analysis

load D140330_Landmark;
% vsTrialPool_E vsTrialPool_L vrDur_E vrDur_L vrDist_E vrDist_L;

% Variables in the landmark struct
% TEOD, EODR, EODA, XH, YH, EODRq, EODRz, ...
% VEL, ACC, ANG, AVEL, ...
% viESAC, vrESAC, vtESAC, vxESAC, vyESAC, img0, ...
% dataID, pathLen_cm, xyFood, xy0, xyStart, ...
% vlZone, Rfood, Rcentre, vrX, vrY, iAnimal, duration
 
%%
get_E = @(x)poolVecFromStruct(vsTrialPool_E, x);
get_L = @(x)poolVecFromStruct(vsTrialPool_L, x);
get_P = @(x)poolVecFromStruct(vsTrialPool_Probe, x);
region_E = @(ml)getRegion(vsTrialPool_E, ml);
region_L = @(ml)getRegion(vsTrialPool_L, ml);
region_P = @(ml)getRegion(vsTrialPool_Probe, ml);

bootCV = @(x)[std(x)/mean(x); bootci(1000, {@(y)std(y)/mean(y), x})];
bootSD = @(x)[std(x); bootci(1000, {@(y)std(y), x})];
bootMean = @(x)[mean(x); bootci(1000, {@(y)mean(y), x})];

%%
rangeLim = [0 5];
threshESAC = 23;
mlMask = getImageMask(img0, rangeLim, 'Food');
strRegion = sprintf('within %d~%d cm of Food', rangeLim(1), rangeLim(2));

figure; 
subplot 221;
vrEODRz_E = get_E('EODRz');
vrEODRz_L = get_L('EODRz');
vrEODRz_P = get_P('EODRz');
plotBarErr(vrEODRz_E(region_E(mlMask)), vrEODRz_L(region_L(mlMask)), vrEODRz_P(region_P(mlMask)));
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title(['EODRz ' strRegion]);

subplot 222;
vrEODA_E = get_E('EODA');
vrEODA_L = get_L('EODA');
vrEODA_P = get_P('EODA');
plotBarErr(vrEODA_E(region_E(mlMask)), ...
           vrEODA_L(region_L(mlMask)), ...
           vrEODA_P(region_P(mlMask)));
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title(['EODA ' strRegion]);

subplot 223;
plotBarErr(vrEODA_E(region_E(mlMask) & vrEODA_E>0) > threshESAC, ...
           vrEODA_L(region_L(mlMask) & vrEODA_L>0) > threshESAC, ...
           vrEODA_P(region_P(mlMask) & vrEODA_P>0) > threshESAC);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title(['ESAC ' strRegion]);

subplot 224;
threshESAC = 23; % 3 SD
plotBarErr(calcDistAsym(vrEODA_E(region_E(mlMask))), ...
           calcDistAsym(vrEODA_L(region_L(mlMask))), ...
           calcDistAsym(vrEODA_P(region_P(mlMask))));
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title(['EODAsym ' strRegion]);
 
 
%% Fig. 4C
dr = .5;
vrR = 0:dr:8;
vrAsym = zeros(size(vrR));
for i=1:numel(vrAsym)
    vrAsym(i) = calcDistAsym(vrEODA_L(region_L(getImageMask(img0, [0 dr] + vrR(i), 'FOOD'))));
end
figure; bar(vrR+dr, vrAsym, 1); axis tight;
xlabel('Dist from food (cm)'); ylabel('EODAsym');
set(gca, 'XTick', 0:2:vrR(end)+dr);


% load mlMaskFood10;
% figure;
% plotBarErr(vrEODRz_E(region_E(mlMask)), vrEODRz_L(region_L(mlMask)));
% set(gca, 'XTickLabel', {'Early', 'Late'});
% title('EODRz within 5-15cm of food');
% 
% mlMask = getImageMask(img0, [5 15], 'FOOD');
% figure; imshow(mlMask);

%% eoda distribution
figure; 
subplot 121; title('Early'); hold on;
for i=1:numel(vsTrialPool_E)
    Y = vsTrialPool_E(i).EODA;
    Y = log10(Y(Y>0));
    ksdensity(Y, 'function', 'survivor');
end
xlabel('log10 EODA pos.');
subplot 122; title('Late'); hold on;
for i=1:numel(vsTrialPool_L)
    Y = vsTrialPool_L(i).EODA;
    Y = log10(Y(Y>0));
    ksdensity(Y, 'function', 'survivor');
end
xlabel('log10 EODA pos.');

figure; 
subplot 121; title('Early'); hold on;
for i=1:numel(vsTrialPool_E)
    Y = vsTrialPool_E(i).EODRz;
%     Y = log10(Y(Y>0));
    ksdensity(Y, 'function', 'survivor');
end
xlabel('EODRz pos.');
subplot 122; title('Late'); hold on;
for i=1:numel(vsTrialPool_L)
    Y = vsTrialPool_L(i).EODRz;
%     Y = log10(Y(Y>0));
    ksdensity(Y, 'function', 'survivor');
end
xlabel('EODRz pos.');

%%
figure; hold on;
for i=1:numel(vsTrialPool_L)
    A = vsTrialPool_L(i).EODA;
    R = vsTrialPool_L(i).EODR;
    DI = -1 .* A(:) ./ R(:).^2;
    ksdensity(DI(DI<0), 'function', 'survivor');
end
title('Late')


%% grid analysis

vrX_E = get_E('vrX');
vrY_E = get_E('vrY');
vrX_L = get_L('vrX');
vrY_L = get_L('vrY');
vrX_P = get_P('vrX');
vrY_P = get_P('vrY');

nGrid = 20;
viX_E = ceil(vrX_E/nGrid);
viY_E = ceil(vrY_E/nGrid);
viX_L = ceil(vrX_L/nGrid);
viY_L = ceil(vrY_L/nGrid);
img0r = imresize(img0a, 1/nGrid);

mrZ_E = nan(size(img0r));
mrZ_L = nan(size(img0r));
vrZ_E = get_E('EODA');
vrZ_L = get_L('EODA');

for iy=1:size(img0r, 1)
    vlY_E = (viY_E == iy);
    vlY_L = (viY_L == iy);
    for ix=1:size(img0r, 2)
        mrZ_E(iy,ix) = mean(vrZ_E(vlY_E & (viX_E == ix)));
        mrZ_L(iy,ix) = mean(vrZ_L(vlY_L & (viX_L == ix)));
    end
end

mlMask = getImageMask(img0, [0 50], 'CENTRE');

figure; 
subplot 121; imshow(rgbmix(img0, imresize(uint8(mrZ_E*255/20), nGrid,'nearest'), mlMask)); title('Early EODA');
subplot 122; imshow(rgbmix(img0, imresize(uint8(mrZ_L*255/20), nGrid,'nearest'), mlMask)); title('Late EODA');

%% back swim density
% [mnVisit1, mnVisit] = calcGridStats(vsTrialPool_E, img0, 'EODAs');
mlMask = getImageMask(img0, [0 50], 'CENTRE');
calcGridStats(vsTrialPool_E, imadjust(img0), 'VEL', @(x)mean(x<0), mlMask); 
title('Early Prob. Backswim ');
figure; imagesc(rand(size(img0))); caxis([0 .5]);


%% EODAs
% [mnVisit1, mnVisit] = calcGridStats(vsTrialPool_E, img0, 'EODAs');
mlMask = getImageMask(img0, [0 50], 'CENTRE');
% calcGridStats(vsTrialPool_Probe, imadjust(img0), 'EODAs', @(x)calcDistAsym(x), mlMask); 
calcGridStats(vsTrialPool_Probe, imadjust(img0), 'EODRz', @(x)mean(x), mlMask); 
title('Probe Mean EODRz');
figure; imagesc(rand(size(img0))); caxis([0 2]);


%% Fig. 2A. Visit density map
pixpercm = 1053.28/(sqrt(2)*100);

region_E = @(ml)getRegion(vsTrialPool_E, ml);
region_L = @(ml)getRegion(vsTrialPool_L, ml);
region_P = @(ml)getRegion(vsTrialPool_P, ml);

img0a = imadjust(img0);

vrImg_E = calcVisitCount(vsTrialPool_E, img0);
vrImg_L = calcVisitCount(vsTrialPool_L, img0);
vrImg_P = calcVisitCount(vsTrialPool_P, img0);

area_E = calcAreaFWHM(vrImg_E) / pixpercm^2; %area in cm^2
area_L = calcAreaFWHM(vrImg_L) / pixpercm^2; %area in cm^2
area_P = calcAreaFWHM(vrImg_P) / pixpercm^2; %area in cm^2

figure;
plotBarErr(area_E, area_L, area_P);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('ESAC prob within 15 cm of food');
ylabel('Search area (1/e of max, 8-conn) cm^2');

%
% vrImg_E1 = zeros(size(img0));
% vrImg_L1 = zeros(size(img0));
% vrImg_E1(sub2ind(size(img0), round(vrY_E(vl_E)), round(vrX_E(vl_E)))) = 1;
% vrImg_L1(sub2ind(size(img0), round(vrY_L(vl_L)), round(vrX_L(vl_L)))) = 1;
% vrImg_E1 = imfilter(vrImg_E1,H);
% vrImg_L1 = imfilter(vrImg_L1,H);
% 
% vrImg_E2 = vrImg_E1 ./ vrImg_E;
% vrImg_L2 = vrImg_L1 ./ vrImg_L;

mlMask = getImageMask(img0, [0 60], 'CENTRE');

figure; 
subplot 221; imshow(rgbmix(img0a, vrImg_E, mlMask)); title('Early Visit counts');
subplot 222; imshow(rgbmix(img0a, vrImg_L, mlMask)); title('Late Visit counts');
subplot 223; imshow(rgbmix(imadjust(img0_P), vrImg_P, mlMask)); title('Probe Visit counts');
% subplot 224; imshow(rgbmix(img0a, vrImg_E2, mlMask)); title('Early Backswim Prob.');
% subplot 225; imshow(rgbmix(img0a, vrImg_L2, mlMask)); title('Late Backswim Prob.');
% subplot 236; imshow(rgbmix(img0a, vrImg_P2, mlMask)); title('Probe Visit counts');
% suptitle('Backswim prob. density density');

%demonstrate FWHM
figure; 
[~, BW1, BW2] = calcAreaFWHM(vrImg_E);
subplot 221; imshow(rgbmix(img0a, vrImg_E, mlMask)); 
subplot 222; imshow(BW1); 
subplot 223; imshow(BW2); 

%% back-swim


vrVEL_E = get_E('VEL');
vrVEL_L = get_L('VEL');
vrVEL_P = get_P('VEL');

figure; 
subplot 221;
mlMask = getImageMask(img0, [0 15], 'LM*');
plotBarErr(vrVEL_E(region_E(mlMask))<0, vrVEL_L(region_L(mlMask))<0, vrVEL_P(region_P(mlMask))<0);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('Prob. Back-swim near landmarks (<15cm)');

subplot 222;
mlMask = getImageMask(img0, [0 15], 'FOOD');
plotBarErr(vrVEL_E(region_E(mlMask))<0, vrVEL_L(region_L(mlMask))<0, vrVEL_P(region_P(mlMask))<0);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('Prob. Back-swim near food (<15cm)');

subplot 223;
mlMask = getImageMask(img0, [0 50], 'CENTRE');
plotBarErr(vrVEL_E(region_E(mlMask))<0, vrVEL_L(region_L(mlMask))<0, vrVEL_P(region_P(mlMask))<0);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('Prob. Back-swim @central zone(<50cm)');

f1 = @(f, m1, m2) sum(f(m1)) / sum(f(m2));
subplot 224;
mlF = getImageMask(img0, [0 15], 'FOOD');
mlA = getImageMask(img0, [0 50], 'CENTRE');
plotBarErr(f1(region_E, mlF, mlA), f1(region_L, mlF, mlA), f1(region_P, mlF, mlA));
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('Fraction of time near food (<15cm)');


%% compare skewness with dsi (distribution symmetry index_

vrSkew_E = [];
vrDSI_E = [];
for i=1:numel(vsTrialPool_E)
    vrSkew_E(end+1) = skewness(poolVecFromStruct(vsTrialPool_E(i), 'EODA'));
    vrDSI_E(end+1) = calcDistAsym(poolVecFromStruct(vsTrialPool_E(i), 'EODA'));
end
vrSkew_L = [];
vrDSI_L = [];
for i=1:numel(vsTrialPool_L)
    vrSkew_L(end+1) = skewness(poolVecFromStruct(vsTrialPool_L(i), 'EODA'));
    vrDSI_L(end+1) = calcDistAsym(poolVecFromStruct(vsTrialPool_L(i), 'EODA'));
end
vrSkew_P = [];
vrDSI_P = [];
for i=1:numel(vsTrialPool_Probe)
    vrSkew_P(end+1) = skewness(poolVecFromStruct(vsTrialPool_Probe(i), 'EODA'));
    vrDSI_P(end+1) = calcDistAsym(poolVecFromStruct(vsTrialPool_Probe(i), 'EODA'));
end

figure; hold on;
plot(vrSkew_E, vrDSI_E, 'b.');
plot(vrSkew_L, vrDSI_L, 'r.');
plot(vrSkew_P, vrDSI_P, 'g.');
xlabel('Skewness');
ylabel('Asymmetry index');


%% Fig. 2A. Location of back-swim

img0a = imadjust(img0);
vl_E = vrVEL_E<0;
vl_L = vrVEL_L<0;
H = fspecial('gaussian', [100 100], 20); %figure; imagesc(H);

vrImg_E = zeros(size(img0));
vrImg_L = zeros(size(img0));
vrImg_E(sub2ind(size(img0), round(vrY_E), round(vrX_E))) = 1;
vrImg_L(sub2ind(size(img0), round(vrY_L), round(vrX_L))) = 1;
vrImg_E = imfilter(vrImg_E,H);
vrImg_L = imfilter(vrImg_L,H);

vrImg_E1 = zeros(size(img0));
vrImg_L1 = zeros(size(img0));
vrImg_E1(sub2ind(size(img0), round(vrY_E(vl_E)), round(vrX_E(vl_E)))) = 1;
vrImg_L1(sub2ind(size(img0), round(vrY_L(vl_L)), round(vrX_L(vl_L)))) = 1;
vrImg_E1 = imfilter(vrImg_E1,H);
vrImg_L1 = imfilter(vrImg_L1,H);

vrImg_E2 = vrImg_E1 ./ vrImg_E;
vrImg_L2 = vrImg_L1 ./ vrImg_L;

load mlMask;
vrImg_E2(~mlMask) = 0;
vrImg_L2(~mlMask) = 0;
% img1 = img0; img1(~mlMask) = 0;

figure; 
subplot 121; imshow(rgbmix(img0a, vrImg_E2, mlMask)); title('Early Backswim Prob.');
subplot 122; imshow(rgbmix(img0a, vrImg_L2, mlMask)); title('Late Backswim Prob.');
suptitle('Backswim prob. density density');


%% fraction of time near food

f1 = @(f, m1, m2) sum(f(m1)) / sum(f(m2));
subplot 224;
mlF = getImageMask(img0, [0 15], 'FOOD');
mlA = getImageMask(img0, [0 50], 'CENTRE');
plotBarErr(f1(region_E, mlF, mlA), f1(region_L, mlF, mlA), f1(region_P, mlF, mlA));
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
title('Fraction of time near food (<15cm)');

