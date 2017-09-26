% Pooled-trial analysis

dname = 'C:\Users\James Jun\Desktop\Landmark_Learning\';

vsFnames = dir([dname '*_Track.mat']);
nFiles = numel(vsFnames);
vcFnames = {};
for i=1:nFiles
    vcFnames{i} = vsFnames(i).name;
end

%% Collect traces
warning off;

vsFishID = {'A', 'B', 'C', 'D'};
for iFish = 1:numel(vsFishID)
    fishID = vsFishID{iFish};
    vsTrial = [];
    for iFile=1:nFiles
        fname = vcFnames{iFile};
        disp(fname);
        if strcmp(fname(4), fishID)
            vsTrial{end+1} = importTrial(load([dname fname]));
        end
    end
    eval(sprintf('vsTrial_%s = vsTrial;', fishID));
end

save D140324_LandmarkGroup vsTrial_*;

% S = makeStruct(TEOD, EODR, EODA, XH, YH, ...
%         VEL, ACC, ANG, AVEL, xy0, ...
%         viESAC, vrESAC, vtESAC, vxESAC, vyESAC, img0, dataID, pathLen_cm);

%%
% load D140324_LandmarkGroup;
load D140330_Landmark;

%% Pool trials
vsFishID = {'A', 'B', 'C', 'D'};
for iFish = 1:numel(vsFishID)
    fishID = vsFishID{iFish};
    eval(sprintf('vsTrial = vsTrial_%s;', fishID));
    vrPathLen = zeros(size(vsTrial));
    vrDuration = zeros(size(vsTrial));
    for iTrial = 1:numel(vsTrial)
        S = vsTrial{iTrial};
        vrPathLen(iTrial) = S.pathLen_cm;
        vrDuration(iTrial) = diff(S.TEOD([1 end]));
    end
    if iFish==1
        img0 = vsTrial{end}.img0; %save the background
    end    
    eval(sprintf('vrPathLen_%s = vrPathLen;', fishID));
    eval(sprintf('vrDuration_%s = vrDuration;', fishID));
end

figure;
subplot 221; bar(log10(vrPathLen_A),1); title('Fish A'); axis([0 60 1.5 4]); ylabel('log10 Dist');
subplot 222; bar(log10(vrPathLen_B),1); title('Fish B'); axis([0 60 1.5 4]); ylabel('log10 Dist');
subplot 223; bar(log10(vrPathLen_C),1); title('Fish C'); axis([0 60 1.5 4]); ylabel('log10 Dist');
subplot 224; bar(log10(vrPathLen_D),1); title('Fish D'); axis([0 60 1.5 4]); ylabel('log10 Dist');

figure;
subplot 221; bar(log10(vrDuration_A),1); title('Fish A'); axis([0 60 0 3]); ylabel('log10 Time');
subplot 222; bar(log10(vrDuration_B),1); title('Fish B'); axis([0 60 0 3]); ylabel('log10 Time');
subplot 223; bar(log10(vrDuration_C),1); title('Fish C'); axis([0 60 0 3]); ylabel('log10 Time');
subplot 224; bar(log10(vrDuration_D),1); title('Fish D'); axis([0 60 0 3]); ylabel('log10 Time');


%% Fig. 1A. Plot Learning curve
% pool 36 groups (4 sessions, 4 individuals = 16 data points)
appendCol = @(x,y)[x, y(:)];
nSessions = 12;

vrPathMu = [];      mrPathCI = [];  vrPathSd = [];
vrDurMu = [];       mrDurCI = [];   vrDurSd = []; 
vrSpeedMu = [];     mrSpeedCI = []; vrSpeedSd = []; 
mrPath = zeros(nSessions, 16);
mrDur = zeros(nSessions, 16);
mrSpeed = zeros(nSessions, 16);
for iSession=1:nSessions
    idxRng = (1:4) + (iSession-1)*4;
    vrPath = [vrPathLen_A(idxRng), vrPathLen_B(idxRng), vrPathLen_C(idxRng), vrPathLen_D(idxRng)];
    vrDur = [vrDuration_A(idxRng), vrDuration_B(idxRng), vrDuration_C(idxRng), vrDuration_D(idxRng)];
    vrSpeed = vrPath ./ vrDur;    
    vrPathMu(end+1) = median(vrPath);
    vrDurMu(end+1) = median(vrDur);
    vrSpeedMu(end+1) = median(vrSpeed);
    vrPathSd(end+1) = std(vrPath);
    vrDurSd(end+1) = std(vrDur);
    vrSpeedSd(end+1) = std(vrSpeed);    
    
    mrPathCI = appendCol(mrPathCI, bootci(1000, {@(x)median(x), vrPath}));
    mrDurCI = appendCol(mrDurCI, bootci(1000, {@(x)median(x), vrDur}));
    mrSpeedCI = appendCol(mrSpeedCI, bootci(1000, {@(x)median(x), vrSpeed}));
    mrPath(iSession, :) = vrPath;
    mrDur(iSession, :) = vrDur;
    mrSpeed(iSession, :) = vrSpeed;
end


figure;
subplot 311; errorbar(1:nSessions, vrPathMu, mrPathCI(1,:), mrPathCI(2,:)); 
    ylabel('Path Length (cm)'); xlabel('Session #'); axis([.5 (nSessions+.5) 0 1500]);
subplot 312; errorbar(1:nSessions, vrDurMu, mrDurCI(1,:), mrDurCI(2,:));  
    ylabel('Duration (s)'); xlabel('Session #'); axis([.5 (nSessions+.5) 0 100]);
subplot 313; errorbar(1:nSessions, vrSpeedMu, mrSpeedCI(1,:), mrSpeedCI(2,:));  
    ylabel('Speed (cm/s)'); xlabel('Session #'); axis([.5 (nSessions+.5) 0 40]);
suptitle('Landmark Learning Group, 4 animals pooled');

figure;
subplot 311; bar(1:nSessions, vrPathSd); 
    ylabel('SD Path Length (cm)'); xlabel('Session #');% axis([.5 (nSessions+.5) 0 1500]);
subplot 312; bar(1:nSessions, vrDurSd);  
    ylabel('SD Duration (s)'); xlabel('Session #'); %axis([.5 (nSessions+.5) 0 100]);
subplot 313; bar(1:nSessions, vrSpeedSd);  
    ylabel('SD Speed (cm/s)'); xlabel('Session #'); %axis([.5 (nSessions+.5) 0 40]);
suptitle('Landmark Learning Group, 4 animals pooled');


%% Fig. 1B. Select trials from early & late learning trials, keep mid 1/8-7/8%
viSession_E = [1:2];
viSession_L = [7:nSessions];
quantLim = [1/8 7/8];
strRange = sprintf('%0.1f~%0.1f%%', quantLim(1)*100, quantLim(2)*100);

figure; 

subplot 311;
[mlPath_E, mlPath_L] = calcHistQuant(log10(mrPath), viSession_E, viSession_L, quantLim, 1.5:.1:4);
xlabel('log10 Distance (cm)');
axis([1.5 4 0 30]); set(gca, {'XTick', 'YTick'}, {1.5:.5:4, 0:10:30});

subplot 312;
% [mlDur_E, mlDur_L] = calcHistQuant(log10(mrDur), viSession_E, viSession_L, quantLim, 0:.15:3.5);
calcHistQuantSelect(log10(mrDur), viSession_E, viSession_L, mlPath_E, mlPath_L, 0:.15:3.5);
xlabel('log10 Duration (s)');
axis([0 3.5 0 30]); set(gca, {'XTick', 'YTick'}, {0:.5:3.5, 0:10:30}); 

subplot 313;
% [mlSpeed_E, mlSpeed_L] = calcHistQuant(mrSpeed, viSession_E, viSession_L, quantLim, [5:1:30]);
calcHistQuantSelect(mrSpeed, viSession_E, viSession_L, mlPath_E, mlPath_L, [5:1:30]);
xlabel('Speed (cm/s)');
axis([5 30 0 20]); set(gca, {'XTick', 'YTick'}, {5:5:30, 0:10:20}); 

figure; 
subplot 121; 
plot(log10(mrDur(:)), log10(mrPath(:)), 'b.');
xlabel('log_{10} Duration (s)'); 
ylabel('log_{10} Distance (cm)');
axis([0 3.5 1.5 4.5]); 
set(gca, 'XTick', 0:1:3.5);
set(gca, 'YTick', 1.5:.5:4.5);

subplot 122; 
plot(log10(mrDur(:)), log10(mrSpeed(:)), 'b.');
xlabel('log_{10} Duration (s)'); 
ylabel('log_{10} Speed (cm/s)');
axis([0 3.5 .75 1.75]); 
set(gca, 'XTick', 0:1:3.5);
set(gca, 'YTick', .75:.25:1.75);

%% Fig. 1D. Early vs. Late bar plot comparison
mlPath_E1 = false(size(mrPath));
mlPath_E1(viSession_E,:) = mlPath_E;
mlPath_L1 = false(size(mrPath));
mlPath_L1(viSession_L,:) = mlPath_L;

bootCV = @(x)[std(x)/mean(x); bootci(1000, {@(y)std(y)/mean(y), x})];
bootSD = @(x)[std(x); bootci(1000, {@(y)std(y), x})];
bootMean = @(x)[mean(x); bootci(1000, {@(y)mean(y), x})];

figure;
subplot 331; plotErrorbar(mrPath, mlPath_E1, mlPath_L1, bootMean, 'Mean Dist. (s)');
subplot 332; plotErrorbar(mrPath, mlPath_E1, mlPath_L1, bootSD, 'SD Dist. (s)');
subplot 333; plotErrorbar(mrPath, mlPath_E1, mlPath_L1, bootCV, 'CV Dist. (s)');
subplot 334; plotErrorbar(mrDur, mlPath_E1, mlPath_L1, bootMean, 'Mean Dur. (s)');
subplot 335; plotErrorbar(mrDur, mlPath_E1, mlPath_L1, bootSD, 'SD Duration (s)');
subplot 336; plotErrorbar(mrDur, mlPath_E1, mlPath_L1, bootCV, 'CV Duration (s)');
subplot 337; plotErrorbar(mrSpeed, mlPath_E1, mlPath_L1, bootMean, 'Mean Speed (s)');
subplot 338; plotErrorbar(mrSpeed, mlPath_E1, mlPath_L1, bootSD, 'SD Speed (s)');
subplot 339; plotErrorbar(mrSpeed, mlPath_E1, mlPath_L1, bootCV, 'CV Speed (s)');


%% FIg. 1E. Mean speed change per animal

% figure; hold on;

mrX = mrSpeed;
vsAnimal = {'A', 'B', 'C', 'D'};
for iAnimal=1:4
    idxRng = 1:4 + 4*(iAnimal-1);
    mlE = false(size(mrSpeed));
    mlL = false(size(mrSpeed));
    mlE(:,idxRng) = mlPath_E1(:,idxRng);
    mlL(:,idxRng) = mlPath_L1(:,idxRng);
    mrEL = [bootMean(mrX(mlE)), bootMean(mrX(mlL))];
    eval(sprintf('mrEL_%s = mrEL;', vsAnimal{iAnimal}));
%     plot(mrEL(1,1), mrEL(1,2), '.');
%     plot(mrEL(1,1)*[1 1], mrEL(2:3,2), '-');
%     plot(mrEL(2:3,1), mrEL(1,2)*[1 1], '-');
end

figure; hold on;
vrX = [(1:4)'-1/8, (1:4)'+1/8];
vrY = [mrEL_A(1,:); mrEL_B(1,:); mrEL_C(1,:); mrEL_D(1,:)];
vrL = [mrEL_A(2,:); mrEL_B(2,:); mrEL_C(2,:); mrEL_D(2,:)];
vrH = [mrEL_A(3,:); mrEL_B(3,:); mrEL_C(3,:); mrEL_D(3,:)];
errorbar(vrX, vrY, vrL, vrH, 'k.');
h = bar(vrY, 1, 'EdgeColor', 'none');
set(h(1), 'FaceColor', 'r');
set(h(2), 'FaceColor', 'b');
set(gca, {'XTick', 'XTickLabel'}, {1:4, {'A', 'B', 'C', 'D'}});
ylabel('Mean swim speed (cm/s)');
set(gca, 'XLim', [.5 4.5]);
set(gca, {'YLim', 'YTick'}, {[0 40], 0:10:40});

%%
% figure;
% qqplot(log(mrPath_E(mlPath_E)), log(mrPath_L(mlPath_L)));
% xlabel('Early'); ylabel('Late');
% 
% figure;
% cdfplot(mrPath_E(:)); hold on;
% cdfplot(mrPath_L(:)); hold on;
% 
% figure;
% cdfplot(mrPath_E(mlPath_E)); hold on;
% cdfplot(mrPath_L(mlPath_L)); hold on;


%% Fig. 2. compare e-saccades distribution between early vs. late learning

appendArray = @(x,y)[x(:); y(:)];
%pool early
vrESAC_E = []; vrEsacEodr_E = []; vrEsacAcc_E = []; vrEsacAvel_E = []; vrEsacVel_E = [];
vrX_E = []; vrY_E = []; vrEODR_E = []; vrEODA_E = []; vrACC_E = []; vrVEL_E = []; vrANG_E = [];
vxESAC_E = []; vyESAC_E = []; vlZone_E = []; vlEsacZone_E = []; vrRfood_E = [];
viAnimal_E = []; viSessionID_E = []; vrEODRq_E = [];
for iSession=viSession_E
    for iAnimal = 1:4
        eval(sprintf('vsTrial = vsTrial_%s;', vsFishID{iAnimal}));    
        for iTrial = 1:4
            if (mlPath_E(iSession-viSession_E(1)+1, iTrial+(iAnimal-1)*4))
                S = vsTrial{iTrial+(iSession-1)*4};
                vlEsacZone_E = appendArray(vlEsacZone_E, S.vlZone(S.viESAC));
                vrESAC_E = appendArray(vrESAC_E, S.vrESAC);
                vrEsacEodr_E = appendArray(vrEsacEodr_E, S.EODR);
                vrEsacAcc_E = appendArray(vrEsacAcc_E, S.ACC);
                vrEsacAvel_E = appendArray(vrEsacAvel_E, S.AVEL);
                vrEsacVel_E = appendArray(vrEsacVel_E, S.VEL);
                vxESAC_E = appendArray(vxESAC_E, S.vrX(S.viESAC));
                vyESAC_E = appendArray(vyESAC_E, S.vrY(S.viESAC));   
                
                vlZone_E = appendArray(vlZone_E, S.vlZone);                
                vrX_E = appendArray(vrX_E, S.vrX);
                vrY_E = appendArray(vrY_E, S.vrY);
                vrEODR_E = appendArray(vrEODR_E, S.EODR);
                vrEODA_E = appendArray(vrEODA_E, S.EODA);
                vrACC_E = appendArray(vrACC_E, S.ACC);
                vrVEL_E = appendArray(vrVEL_E, S.VEL);
                vrEODRq_E = appendArray(vrEODRq_E, S.EODRq);
                
                viSessionID_E = appendArray(viSessionID_E, iSession * ones(size(S.EODR)));
                viAnimal_E = appendArray(viAnimal_E, iAnimal * ones(size(S.EODR)));
            end
        end
    end
end

%pool late
vrESAC_L = []; vrEsacEodr_L = []; vrEsacAcc_L = []; vrEsacAvel_L = []; vrEsacVel_L = [];
vrX_L = []; vrY_L = []; vrEODR_L = []; vrEODA_L = []; vrACC_L = []; vrVEL_L = []; vrANG_L = [];
vxESAC_L = []; vyESAC_L = []; vlZone_L = []; vlEsacZone_L = []; vrRfood_L = [];
viAnimal_L = []; viSessionID_L = []; vrEODRq_L = [];
for iSession=viSession_L
    for iAnimal = 1:4
        eval(sprintf('vsTrial = vsTrial_%s;', vsFishID{iAnimal}));    
        for iTrial = 1:4
            if (mlPath_L(iSession-viSession_L(1)+1, iTrial+(iAnimal-1)*4))
                S = vsTrial{iTrial+(iSession-1)*4};
                vlEsacZone_L = appendArray(vlEsacZone_L, S.vlZone(S.viESAC));
                vrESAC_L = appendArray(vrESAC_L, S.vrESAC);
                vrEsacEodr_L = appendArray(vrEsacEodr_L, S.EODR);
                vrEsacAcc_L = appendArray(vrEsacAcc_L, S.ACC);
                vrEsacAvel_L = appendArray(vrEsacAvel_L, S.AVEL);
                vrEsacVel_L = appendArray(vrEsacVel_L, S.VEL);
                vxESAC_L = appendArray(vxESAC_L, S.vrX(S.viESAC));
                vyESAC_L = appendArray(vyESAC_L, S.vrY(S.viESAC));
                
                vlZone_L = appendArray(vlZone_L, S.vlZone);                
                vrX_L = appendArray(vrX_L, S.vrX);
                vrY_L = appendArray(vrY_L, S.vrY);
                vrEODR_L = appendArray(vrEODR_L, S.EODR);
                vrEODA_L = appendArray(vrEODA_L, S.EODA);
                vrACC_L = appendArray(vrACC_L, S.ACC);
                vrVEL_L = appendArray(vrVEL_L, S.VEL);
                vrRfood_L = appendArray(vrRfood_L, S.Rfood);
                vrEODRq_L = appendArray(vrEODRq_L, S.EODRq);
                
                viSessionID_L = appendArray(viSessionID_L, iSession * ones(size(S.EODR)));
                viAnimal_L = appendArray(viAnimal_L, iAnimal * ones(size(S.EODR)));                
            end
        end
    end
end
vlEsacZone_E = logical(vlEsacZone_E);
vlEsacZone_L = logical(vlEsacZone_L);
vlZone_E = logical(vlZone_E);
vlZone_L = logical(vlZone_L);
viX_E = round(vrX_E); viY_E = round(vrY_E);
viX_L = round(vrX_L); viY_L = round(vrY_L);

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
subplot 321; imshow(rgbmix(img0a, vrImg_E, mlMask)); title('Early Visit counts');
subplot 322; imshow(rgbmix(img0a, vrImg_L, mlMask)); title('Late Visit counts');
subplot 323; imshow(rgbmix(img0a, vrImg_E1, mlMask)); title('Early Backswim counts');
subplot 324; imshow(rgbmix(img0a, vrImg_L1, mlMask)); title('Late Backswim counts');
subplot 325; imshow(rgbmix(img0a, vrImg_E2, mlMask)); title('Early Backswim Prob.');
subplot 326; imshow(rgbmix(img0a, vrImg_L2, mlMask)); title('Late Backswim Prob.');
suptitle('Backswim prob. density density');


%% Fig. 2B. Backswim prob. near landmarks and food location. bar plot comparison
getMatVal = @(mlMask, viY, viX) mlMask(sub2ind(size(mlMask), viY, viX));
% load mlMaskLandmark10;
load mlMaskFood20;
vlVisit_E = getMatVal(mlMask, viY_E, viX_E);
vlVisit_L = getMatVal(mlMask, viY_L, viX_L);
vlBack_E = vrVEL_E(vlVisit_E) < 0;
vlBack_L = vrVEL_L(vlVisit_L) < 0;
vpBack_LM_E = bootMean(vlBack_E);
vpBack_LM_L = bootMean(vlBack_L);
mr_LM = [vpBack_LM_E, vpBack_LM_L];

figure; 
bar([1 2], mr_LM(1,:), .5);
hold on; 
errorbar([1 2], mr_LM(1,:), mr_LM(2,:), mr_LM(3,:), 'r.');
ylabel('Prob. Backswim');
set(gca, 'XTickLabel', {'Early', 'Late'})


%% Fig. 2C. Backswim prob. change. todo: per each animal


% backswim onset detection
calcProbBack = @(y)[mean(y < 0); bootci(1000, {@(x)mean(x < 0), y}, 'type', 'cper')];
vec_E = calcProbBack(vrVEL_E(vlZone_E));
vec_L = calcProbBack(vrVEL_L(vlZone_L));

figure; 
bar([1 2], [vec_E(1), vec_L(1)], .5);
hold on; 
errorbar([1 2], [vec_E(1), vec_L(1)], [vec_E(2), vec_L(2)], [vec_E(3), vec_L(3)], 'r.');
ylabel('Prob. Backswim');
set(gca, 'XTickLabel', {'Early', 'Late'})
% [h p] = kstest2(calcBootEsac(vrESAC_E), calcBootEsac(vrESAC_L))

%% Fig. 2E. Visit density near food given region of interest. CDF plot
bootMean = @(x)[mean(x); bootci(100, {@(y)mean(y), x}, 'Type', 'cper')];

load mlMaskCentre50; %active zone location
calcDistCM = @(X, Y, xy0)sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;

viX_E = round(vrX_E); viY_E = round(vrY_E);
viX_L = round(vrX_L); viY_L = round(vrY_L);
vlVisit_E = getMatVal(mlMask, viY_E, viX_E);
vlVisit_L = getMatVal(mlMask, viY_L, viX_L);
vrR_E = vrRfood_E(vlVisit_E);
vrR_L = vrRfood_L(vlVisit_L);

vrDist = 0:2.5:60;
mpVisit_E = [];
mpVisit_L = [];
for i=1:numel(vrDist)
    mpVisit_E = [mpVisit_E, bootMean(vrR_E > vrDist(i))];
    mpVisit_L = [mpVisit_L, bootMean(vrR_L > vrDist(i))];
end
figure; plot(vrDist, mpVisit_E, 'r');
hold on; plot(vrDist, mpVisit_L, 'b');

%% Fig. 2D. bar plot food dist
% Bar plot
vec_E = bootMean(vrR_E <10);
vec_L = bootMean(vrR_L <10);
figure; 
bar([1 2], [vec_E(1), vec_L(1)], .5);
hold on; 
errorbar([1 2], [vec_E(1), vec_L(1)], [vec_E(2), vec_L(2)], [vec_E(3), vec_L(3)], 'r.');
ylabel('Fraction near food (within 10 cm)');
set(gca, 'XTickLabel', {'Early', 'Late'})

%% -----------------------------------------------------------------------
%% Fig. 3E. EODA location map

thresh = std([vrEODA_E; vrEODA_L])*2;

img0a = imadjust(img0);
H = fspecial('gaussian', [100 100], 20); %figure; imagesc(H);
vrImg_E = zeros(size(img0));
vrImg_L = zeros(size(img0));
vrImg_E(sub2ind(size(img0), round(vrY_E(vrEODA_E>0)), round(vrX_E(vrEODA_E>0)))) = 1;
vrImg_L(sub2ind(size(img0), round(vrY_L(vrEODA_L>0)), round(vrX_L(vrEODA_L>0)))) = 1;
vrImg_E = imfilter(vrImg_E,H);
vrImg_L = imfilter(vrImg_L,H);

vrImg_E1 = zeros(size(img0));
vrImg_L1 = zeros(size(img0));
vrImg_E1(sub2ind(size(img0), round(vrY_E(vrEODA_E>thresh)), round(vrX_E(vrEODA_E>thresh)))) = 1;
vrImg_L1(sub2ind(size(img0), round(vrY_L(vrEODA_L>thresh)), round(vrX_L(vrEODA_L>thresh)))) = 1;
vrImg_E1 = imfilter(vrImg_E1,H);
vrImg_L1 = imfilter(vrImg_L1,H);

vrImg_E2 = vrImg_E1 ./ vrImg_E;
vrImg_L2 = vrImg_L1 ./ vrImg_L;

load mlMask;
vrImg_E2(~mlMask) = 0;
vrImg_L2(~mlMask) = 0;
% img1 = img0; img1(~mlMask) = 0;

figure; 
% subplot 321; imshow(rgbmix(img0a, vrImg_E, mlMask)); title('Early Visit counts');
% subplot 322; imshow(rgbmix(img0a, vrImg_L, mlMask)); title('Late Visit counts');
subplot 221; imshow(rgbmix(img0a, vrImg_E1, mlMask)); title('Early ESAC map');
subplot 222; imshow(rgbmix(img0a, vrImg_L1, mlMask)); title('Late ESAC map');
subplot 223; imshow(rgbmix(img0a, vrImg_E2, mlMask)); title('Early ESAC prob.');
subplot 224; imshow(rgbmix(img0a, vrImg_L2, mlMask)); title('Late ESAC prob.');
suptitle('EODA outlier prob.');


%% Fig. 3F. ESAC prob. 

thresh = std([vrEODA_E; vrEODA_L])*2;

getMatVal = @(mlMask, viY, viX) mlMask(sub2ind(size(mlMask), viY, viX));
load mlMaskFood10;
% load mlMaskLandmark10;
vlVisit_E = getMatVal(mlMask, viY_E, viX_E);
vlVisit_L = getMatVal(mlMask, viY_L, viX_L);
vlBack_E = vrEODA_E(vlVisit_E & vrEODA_E > 0) > thresh;
vlBack_L = vrEODA_E(vlVisit_L & vrEODA_L > 0) > thresh;
vpBack_LM_E = bootMean(vlBack_E);
vpBack_LM_L = bootMean(vlBack_L);
mr_LM = [vpBack_LM_E, vpBack_LM_L];

figure; 
bar([1 2], mr_LM(1,:), .5);
hold on; 
errorbar([1 2], mr_LM(1,:), mr_LM(2,:), mr_LM(3,:), 'r.');
ylabel('Prob. ESAC');
set(gca, 'XTickLabel', {'Early', 'Late'})
title('Prob. ESAC near Food');

%% Fig. 3B. qq plot for EODA
limInset = [1.4 1.8];
figure; 
subplot 121; 
qqplot(log10(vrEODA_E(vlZone_E & vrEODA_E>0)), log10(vrEODA_L(vlZone_L & vrEODA_L>0)));
xlabel('log10 EODA Early (Hz/s)'); ylabel('log10 EODA Late (Hz/s)');
axis([-4 3 -4 3]); axis square;

subplot 122;
qqplot(log10(vrEODA_E(vlZone_E & vrEODA_E>0)), log10(vrEODA_L(vlZone_L & vrEODA_L>0)));
xlabel('log10 EODA Early (Hz/s)'); ylabel('log10 EODA Late (Hz/s)');
axis([limInset, limInset]); axis square; 
set(gca, {'XTick', 'YTick'}, {limInset(1):.1:limInset(2), limInset(1):.1:limInset(2)});


% set(gca, {'XTick', 'YTick'}, {0:20:100, 0:20:100}); axis([0 100 0 100]); axis square;

% subplot 122; qqplot(log10(vrEODA_E(vlZone_E & vrEODA_E>0)), log10(vrEODA_L(vlZone_L & vrEODA_L>0)));
% xlabel('EODA Early (Hz/s)'); ylabel('EODA Late (Hz/s)');
% set(gca, {'XScale', 'YScale'}, {'log', 'log'})
% set(gca, {'XTick', 'YTick'}, {limInset(1):5:limInset(2), limInset(1):5:limInset(2)}); axis([limInset, limInset]); axis square;

figure; plotKS2(log10(vrEODA_E(vlZone_E & vrEODA_E>0)), log10(vrEODA_L(vlZone_L & vrEODA_L>0)));
% figure; plotKS2(vrESAC_E(vlEsacZone_E), vrESAC_L(vlEsacZone_L));
% quantile(vrESAC_E(vlEsacZone_E), .92)

%% Fig. 3C. EODA distribution (positive only)

figure; hold on;
[xi, vrF, vrF_L, vrF_H] = plotDistBootci(vrEODA_E(vlZone_E & vrEODA_E>0));
plot(xi, [vrF, vrF_L, vrF_H], 'b');
[xi, vrF, vrF_L, vrF_H] = plotDistBootci(vrEODA_L(vlZone_L & vrEODA_L>0));
plot(xi, [vrF, vrF_L, vrF_H], 'r');

% figure; ksdensity(vrESAC_E, 'function', 'survivor'); 
% hold on; ksdensity(vrESAC_L, 'function', 'survivor'); 
legend({'Early', 'Late'});
xlabel('EODA (Hz/s)');
ylabel('Survivor Function');
% set(gca, 'YTick', 0:.01:.05);
set(gca, 'YScale', 'log'); 
set(gca, 'XScale', 'linear');
axis([0 50 1e-4 1]);
% plot(threshESAC*[1 1], get(gca, 'YLim'), 'k');

% figure; cdfplot(log10(vrESAC_L)); hold on; cdfplot(log10(vrESAC_E)); xlabel('log10 ESAC');
% legend({'Late', 'Early'});


%% Fig. 3D. EODA dist. during backswim
thresh = std([vrEODA_E(vrEODA_E>0); vrEODA_L(vrEODA_L>0)])*2;

load mlMaskCentre50Food5;
vlLoc_E = getMatVal(mlMask, viY_E, viX_E);
vlLoc_L = getMatVal(mlMask, viY_L, viX_L);
vlESAC_E = vrEODA_E(vlLoc_E & vrEODA_E>0) > thresh;
vlESAC_L = vrEODA_L(vlLoc_L & vrEODA_L>0) > thresh;
vpE = bootMean(vlESAC_E);
vpL = bootMean(vlESAC_L);
mr = [vpE, vpL];

figure; 
bar([1:2], mr(1,:), .5);
hold on; 
errorbar([1:2], mr(1,:), mr(2,:), mr(3,:), 'r.');
ylabel(sprintf('Prob. ESAC', threshEODA));
set(gca, 'XTickLabel', {'Early', 'Late'})

figure; hold on;
[xi, vrF, vrF_L, vrF_H] = plotDistBootci(vrEODA_E(vlLoc_E & vrEODA_E>0));
plot(xi, [vrF, vrF_L, vrF_H], 'r');
[xi, vrF, vrF_L, vrF_H] = plotDistBootci(vrEODA_L(vlLoc_L & vrEODA_L>0));
plot(xi, [vrF, vrF_L, vrF_H], 'b');
legend({'Early', 'Late'});
xlabel('EODA (Hz/s)');
ylabel('Survivor Function');
% set(gca, 'YTick', 0:.01:.05);
set(gca, 'YScale', 'log'); 
set(gca, 'XScale', 'linear');
axis([0 50 1e-4 1]);


%%----------------------------------------------------------------------
%% Fig. 4A. EODR distribution change

%quantile-normalize. create a quantile vector vqEODR_L and vqEODR_E
%normalized by the same session

mrEODR_med_L = zeros(4, 12);
for iAnimal = 1:4    
    for iSession=viSession_L
        vr = vrEODRq_L(viAnimal_L == iAnimal & viSessionID_L == iSession);
        mrEODR_med_L(iAnimal, iSession) = median(vr);
    end
end

figure; imagesc(mrEODR_med_L(1:4, 7:12));

%% Query location ofo quantile shift

img0a = imadjust(img0);
H = fspecial('gaussian', [100 100], 20); %figure; imagesc(H);
vrImg_E = zeros(size(img0));
vrImg_L = zeros(size(img0));
vrImg_E(sub2ind(size(img0), viY_E, viX_E)) = 1;
vrImg_L(sub2ind(size(img0), viY_L, viX_L)) = 1;
vrImg_E = imfilter(vrImg_E,H);
vrImg_L = imfilter(vrImg_L,H);

vrImg_E1 = zeros(size(img0));
vrImg_L1 = zeros(size(img0));
vrImg_E1(sub2ind(size(img0), viY_E(vrEODRq_E>0), viX_E(vrEODRq_E>0))) = 1;
vrImg_L1(sub2ind(size(img0), viY_L(vrEODRq_L>0), viX_L(vrEODRq_L>0))) = 1;
vrImg_E1 = imfilter(vrImg_E1,H);
vrImg_L1 = imfilter(vrImg_L1,H);

vrImg_E2 = vrImg_E1;% ./ vrImg_E;
vrImg_L2 = vrImg_L1;% ./ vrImg_L;

load mlMask;
vrImg_E2(~mlMask) = 0;
vrImg_L2(~mlMask) = 0;
% img1 = img0; img1(~mlMask) = 0;

figure; 
% subplot 321; imshow(rgbmix(img0a, vrImg_E, mlMask)); title('Early Visit counts');
% subplot 322; imshow(rgbmix(img0a, vrImg_L, mlMask)); title('Late Visit counts');
subplot 121; imshow(rgbmix(img0a, uint8(vrImg_E2*255), mlMask)); title('Early prob. EODR Quant > .5');
subplot 122; imshow(rgbmix(img0a, uint8(vrImg_L2*255), mlMask)); title('Late prob. EODR Quant > .5');
suptitle('EODA outlier prob.');

%% quantile shift location
% vrEODRq_E1 = calcQuantileShift(vrEODR_E);
% vrEODRq_L1 = calcQuantileShift(vrEODR_L);
subtMu = @(x)x-mean(x);
subtAbsMu = @(x)abs(x)-mean(abs(x));

load mlMaskFood10;
% load mlMaskLandmark10;
% load mlMaskCentre50;
vlLoc_E = getMatVal(mlMask, viY_E, viX_E);
vlLoc_L = getMatVal(mlMask, viY_L, viX_L);
vrEODRq_E1 = vrEODRq_E(vlLoc_E);
vrEODRq_L1 = vrEODRq_L(vlLoc_L);
vrVEL_E1 = vrVEL_E(vlLoc_E);
vrVEL_L1 = vrVEL_L(vlLoc_L);

ncorr = 1000;
vrTplot = (-ncorr:ncorr)/100;
figure; hold on; title('EODA vs Accel');
plot(vrTplot, xcorr(subtMu(vrEODA_E(vlLoc_E)),subtMu(vrACC_E(vlLoc_E)), ncorr, 'coeff'),'r'); xlabel('Time of EODA (s)');
plot(vrTplot, xcorr(subtMu(vrEODA_L(vlLoc_L)),subtMu(vrACC_L(vlLoc_L)), ncorr, 'coeff'),'b'); xlabel('Time of EODA (s)');
% set(gca, 'XScale', 'log');

figure; plot(-1000:1000, xcorr(subtMu(vrACC_E(vlLoc_E)),subtMu(vrEODA_E(vlLoc_E)), 1000, 'coeff')); xlabel('Time of ACC');

vlEODRq_E = vrEODRq_E(vlLoc_E)>0;
vlEODRq_L = vrEODRq_L(vlLoc_L)>0;
vpE = bootMean(vlEODRq_E);
vpL = bootMean(vlEODRq_L);
mr = [vpE, vpL];

figure; 
bar([1:2], mr(1,:), .5);
hold on; 
errorbar([1:2], mr(1,:), mr(2,:), mr(3,:), 'r.');
ylabel(sprintf('Mean EODR Quantile shift', threshEODA));
set(gca, 'XTickLabel', {'Early', 'Late'})