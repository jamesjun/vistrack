% Pooled-trial analysis

dname = 'C:\Users\James Jun\Desktop\New folder\';

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
load D140324_LandmarkGroup;

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

nSessions = 12;

vrPathMu = [];  vrPath_L = [];  vrPath_H = [];
vrDurMu = [];   vrDur_L = [];   vrDur_H = [];
mrPath = zeros(nSessions, 16);
mrDur = zeros(nSessions, 16);
for iSession=1:nSessions
    idxRng = (1:4) + (iSession-1)*4;
    vrPath = [vrPathLen_A(idxRng), vrPathLen_B(idxRng), vrPathLen_C(idxRng), vrPathLen_D(idxRng)];
    vrDur = [vrDuration_A(idxRng), vrDuration_B(idxRng), vrDuration_C(idxRng), vrDuration_D(idxRng)];
    vrPathMu(end+1) = median(vrPath);
    vrDurMu(end+1) = median(vrDur);
    a = bootci(1000, {@(x)median(x), vrPath});
    vrPath_L(end+1) = a(1);
    vrPath_H(end+1) = a(2);
    a = bootci(1000, {@(x)median(x), vrDur});
    vrDur_L(end+1) = a(1);
    vrDur_H(end+1) = a(2);
    mrPath(iSession, :) = vrPath;
    mrDur(iSession, :) = vrDur;
end


figure;
subplot 211; errorbar(1:nSessions, vrPathMu, vrPath_L, vrPath_H); 
    ylabel('Path Length (cm)'); xlabel('Session #'); axis([.5 (nSessions+.5) 0 1500]);
subplot 212; errorbar(1:nSessions, vrDurMu, vrDur_L, vrDur_H); 
    ylabel('Duration (s)'); xlabel('Session #'); axis([.5 (nSessions+.5) 0 100]);
suptitle('Landmark Learning Group, 4 animals pooled');


%% Fig. 1B. Select trials from early & late learning trials, keep mid 25-75%
viSession_E = [1:3];
viSession_L = [4:nSessions];
quantLim = [.25 .75];
strRange = sprintf('%d~%d%%', round(quantLim(1)*100), round(quantLim(2)*100));

mrPath_E = mrPath(viSession_E,:);
mrPath_L = mrPath(viSession_L,:);
pathQuant_E = quantile(mrPath_E(:), quantLim);
pathQuant_L = quantile(mrPath_L(:), quantLim);
mlPath_E = (mrPath_E >= pathQuant_E(1)) & (mrPath_E < pathQuant_E(2));
mlPath_L = (mrPath_L >= pathQuant_L(1)) & (mrPath_L < pathQuant_L(2));
pathLenTot_E = sum(sum(mrPath_E(mlPath_E)));
pathLenTot_L = sum(sum(mrPath_L(mlPath_L)));

vrLogPath_E = log10(mrPath_E(:));
vrLogPath_L = log10(mrPath_L(:));
xi = 1.5:.1:4;
vnE = hist(vrLogPath_E, xi); 
vnL = hist(vrLogPath_L, xi);
vnE1 = hist(vrLogPath_E(mlPath_E), xi);
vnL1 = hist(vrLogPath_L(mlPath_L), xi);

figure; 
subplot 121;    bar(xi, [vnE; vnL]', 1, 'EdgeColor', 'none');
legend({'Early all', 'Late all'});
xlabel('log10 PathLen (cm)');   ylabel('# Obs.');

subplot 122;    bar(xi, [vnE1; vnL1]', 1, 'EdgeColor', 'none');
legend({['Early ' strRange], ['Late ' strRange]});
xlabel('log10 PathLen (cm)');   ylabel('# Obs.');

figure;
qqplot(log(mrPath_E(mlPath_E)), log(mrPath_L(mlPath_L)));
xlabel('Early'); ylabel('Late');

figure;
cdfplot(mrPath_E(:)); hold on;
cdfplot(mrPath_L(:)); hold on;

figure;
cdfplot(mrPath_E(mlPath_E)); hold on;
cdfplot(mrPath_L(mlPath_L)); hold on;


% compare e-saccades distribution between early vs. late learning

appendArray = @(x,y)[x(:); y(:)];
%pool early
vrESAC_E = []; vrEsacEodr_E = []; vrEsacAcc_E = []; vrEsacAvel_E = []; vrEsacVel_E = [];
vrX_E = []; vrY_E = []; vrEODR_E = []; vrEODA_E = []; vrACC_E = []; vrVEL_E = []; vrANG_E = [];
vxESAC_E = []; vyESAC_E = []; vlZone_E = []; vlEsacZone_E = []; vrRfood_E = [];
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
                vrRfood_E = appendArray(vrRfood_E, S.Rfood);
            end
        end
    end
end

%pool late
vrESAC_L = []; vrEsacEodr_L = []; vrEsacAcc_L = []; vrEsacAvel_L = []; vrEsacVel_L = [];
vrX_L = []; vrY_L = []; vrEODR_L = []; vrEODA_L = []; vrACC_L = []; vrVEL_L = []; vrANG_L = [];
vxESAC_L = []; vyESAC_L = []; vlZone_L = []; vlEsacZone_L = []; vrRfood_L = [];
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
            end
        end
    end
end
vlEsacZone_E = logical(vlEsacZone_E);
vlEsacZone_L = logical(vlEsacZone_L);
vlZone_E = logical(vlZone_E);
vlZone_L = logical(vlZone_L);

%% Fig. 1D. qq plot for EODA
limInset = [25 40];
figure; 
subplot 121; qqplot(vrEODA_E(vlZone_E & vrEODA_E>0), vrEODA_L(vlZone_L & vrEODA_L>0));
xlabel('EODA Early (Hz/s)'); ylabel('EODA Late (Hz/s)');
set(gca, {'XTick', 'YTick'}, {0:20:100, 0:20:100}); axis([0 100 0 100]); axis square;

subplot 122; qqplot(vrEODA_E(vlZone_E & vrEODA_E>0), vrEODA_L(vlZone_L & vrEODA_L>0));
xlabel('EODA Early (Hz/s)'); ylabel('EODA Late (Hz/s)');
set(gca, {'XTick', 'YTick'}, {limInset(1):5:limInset(2), limInset(1):5:limInset(2)}); axis([limInset, limInset]); axis square;

% figure; plotKS2(vrESAC_E(vlEsacZone_E), vrESAC_L(vlEsacZone_L));
% quantile(vrESAC_E(vlEsacZone_E), .92)

%% Fig. 1E. EODA distribution (positive only)

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
axis([0 60 1e-4 1]);
set(gca, 'XLim', [0 50]);
% plot(threshESAC*[1 1], get(gca, 'YLim'), 'k');

% figure; cdfplot(log10(vrESAC_L)); hold on; cdfplot(log10(vrESAC_E)); xlabel('log10 ESAC');
% legend({'Late', 'Early'});


%% Fig. 1F,G,H. EODA dist. during backswim
threshEODA = 35;

% backswim onset detection
calcProbOL = @(y)[mean(y>threshEODA); bootci(1000, {@(x)mean(x>threshEODA), y}, 'type', 'cper')];

vec_E = calcProbOL(vrEODA_E(vlZone_E & vrEODA_E>0));
vec_L = calcProbOL(vrEODA_L(vlZone_L & vrEODA_L>0));
mr = [vec_E(:), vec_L(:)];
figure; 
bar([1:2], mr(1,:), .5);
hold on; 
errorbar([1:2], mr(1,:), mr(2,:), mr(3,:), 'r.');
ylabel(sprintf('Prob. EODA > %0.0f Hz/s', threshEODA));
set(gca, 'XTickLabel', {'Early', 'Late'})


vecF_E = calcProbOL(vrEODA_E(vrVEL_E>0 & vlZone_E & vrEODA_E>0));
vecB_E = calcProbOL(vrEODA_E(vrVEL_E<0 & vlZone_E & vrEODA_E>0));
vecF_L = calcProbOL(vrEODA_L(vrVEL_L>0 & vlZone_L & vrEODA_L>0));
vecB_L = calcProbOL(vrEODA_L(vrVEL_L<0 & vlZone_L & vrEODA_L>0));
mr = [vecF_E(:), vecB_E(:), vecF_L(:), vecB_L(:)];
figure; 
bar([1:4], mr(1,:), .5);
hold on; 
errorbar([1:4], mr(1,:), mr(2,:), mr(3,:), 'r.');
ylabel(sprintf('Prob. EODA > %0.0f Hz/s', threshEODA));
set(gca, 'XTickLabel', {'Early-F', 'Early-B', 'Late-F', 'Late-B'})
% [h p] = kstest2(calcBootEsac(vrESAC_E), calcBootEsac(vrESAC_L))


%--------------------------------------------------------------------------
%% Fig. 2A. position-dependent stats
calcProbOL = @(y)mean(y>20);
    
vrProbOL_E = [];
vrProbOL_L = [];
vrDist = 5:2:15;
for i=1:numel(vrDist)-1
    vrProbOL_E(end+1) = calcProbOL(vrEODA_E(vrEODA_E>0 & vrRfood_E >= vrDist(i) & vrRfood_E < vrDist(i+1)));
    vrProbOL_L(end+1) = calcProbOL(vrEODA_L(vrEODA_L>0 & vrRfood_L >= vrDist(i) & vrRfood_L < vrDist(i+1)));
end
vrX = vrDist(1:end-1) + 2/2;
figure;  AX = [];
subplot 121; bar(vrX, vrProbOL_E); title('Early'); AX(1) = gca;
subplot 122; bar(vrX, vrProbOL_L); title('Late'); AX(2) = gca;
linkaxes(AX, 'xy');

%% show correlation
figure;
X = vrVEL_L(vlZone_L & vrEODA_L>0 & vrVEL_L>0);
Y = vrEODA_L(vlZone_L & vrEODA_L>0 & vrVEL_L>0);
figure; plot(X, Y, '.'); 
xlabel('Vel.'); ylabel('EODA');


%% Triggered average of EODR before backswim

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


%% query top 5% quantile
quantile(abs(vrEODA_E(vlZone_E)), .99)
quantile(abs(vrEODA_L(vlZone_L)), .99)

quantile(vrEODA_E(vlZone_E & vrEODA_E>0), .99)
quantile(vrEODA_L(vlZone_L & vrEODA_L>0), .99)
figure; qqplot(vrEODA_E(vlZone_E & vrEODA_E>0), vrEODA_L(vlZone_L & vrEODA_L>0));

figure; ksdensity(vrEODA_E(vlZone_E & vrEODA_E>0), 'function', 'pdf');
hold on; ksdensity(vrEODA_L(vlZone_L & vrEODA_L>0), 'function', 'pdf');
legend({'Early', 'Late'});
axis([0 70 1e-5 1e0]);


quantile(vrESAC_E(vlEsacZone_E), .99);
quantile(vrESAC_L(vlEsacZone_L), .99);

threshESAC = 20; %2SDgreater
calcEsacCnt = @(y)[mean(((y >= threshESAC))); bootci(1000, {@(x)mean((x >= threshESAC)), y})];
% bootEsacCnt = @(y)bootstrp(1000, @(x)sum(x >= threshESAC), y);
vec_E = calcEsacCnt(vrESAC_E(vlEsacZone_E));% / sum(vlZone_E)/100;
vec_L = calcEsacCnt(vrESAC_L(vlEsacZone_L));% / sum(vlZone_L)/100;

% vrESAC95_E = vrESAC_E(vrESAC_E >= quantile(vrESAC_E, .95));
% vrESAC95_L = vrESAC_L(vrESAC_L >= quantile(vrESAC_L, .95));
% calcMed = @(y)[median(y); bootci(1000, {@(x)median(x), y})];
% vec_E = calcMed(vrESAC_E);
% vec_L = calcMed(vrESAC_L);
figure; 
bar([1 2], [vec_E(1), vec_L(1)], .5);
hold on; 
errorbar([1 2], [vec_E(1), vec_L(1)], [vec_E(2), vec_L(2)], [vec_E(3), vec_L(3)], 'r.');
ylabel(sprintf('# ESAC (>%0.1f) per sec. (#/s)', threshESAC));
set(gca, 'XTickLabel', {'Early', 'Late'})
% [h p] = kstest2(calcBootEsac(vrESAC_E(vlEsacZone_E)), calcBootEsac(vrESAC_L(vlEsacZone_L)))



%% Plot the acceleration-triggered E-Sac

figure;
subplot 221; plotBox(vrESAC_E(vlEsacZone_E), abs(vrEsacAcc_E(vlEsacZone_E)), [0 40 1000]); xlabel('ESAC'); ylabel('|Accel|'); title('Accel vs. E-Sac. Early');
subplot 222; plotBox(vrESAC_E(vlEsacZone_E), abs(vrEsacVel_E(vlEsacZone_E)), [0 40 1000]); xlabel('ESAC'); ylabel('abs Vel.'); title('Accel vs. E-Sac. Early');
subplot 223; plotBox(vrESAC_L(vlEsacZone_L), abs(vrEsacAcc_L(vlEsacZone_L)), [0 40 1000]); xlabel('ESAC'); ylabel('|Accel|'); title('Accel vs. E-Sac. Late');
subplot 224; plotBox(vrESAC_L(vlEsacZone_L), abs(vrEsacVel_L(vlEsacZone_L)), [0 40 1000]); xlabel('ESAC'); ylabel('abs Vel.'); title('Accel vs. E-Sac. Late');

%% locations associated with ESAC: early vs. late. rotate and pool on the picture. 

figure;
subplot 121; plotBox(double(vrVEL_E(vlZone_E)>0), abs(vrEODA_E(vlZone_E)), [-1 .5 1]); 
ylabel('|EODA| (Hz/s)'); title('Early');  grid on; set(gca, 'XTick', [1 2]);
axis([.5 2.5 -20 20]);
set(gca, 'XTickLabel', {'Back', 'Forward'}); xlabel('Swim direction');

subplot 122; plotBox(double(vrVEL_L(vlZone_L)>0), abs(vrEODA_L(vlZone_L)), [-1 .5 1]); 
ylabel('|EODA| (Hz/s)'); title('Late'); grid on; set(gca, 'XTick', [1 2]);
axis([.5 2.5 -20 20]); 
set(gca, 'XTickLabel', {'Back', 'Forward'}); xlabel('Swim direction');

% [h p] = ttest2(vrEODR_L(vrVEL_L>0), vrEODR_L(vrVEL_L<0))

%% region-based analysis. location of e-saccades
threshESAC = 15;

vlESAC = vrESAC_E > threshESAC;
vxESAC1_E = vxESAC_E(vlESAC );
vyESAC1_E = vyESAC_E(vlESAC );

vlESAC = vrESAC_L > threshESAC;
vxESAC1_L = vxESAC_L(vlESAC );
vyESAC1_L = vyESAC_L(vlESAC );

figure; imshow(img0); hold on; title('Loc. E-Sacc > 40 during Early learning');
plot(vxESAC1_E, vyESAC1_E, 'b.');
plot(vxESAC1_L, vyESAC1_L, 'r.');

% e-saccades spread function
vrImg_E = zeros(size(img0));
vrImg_L = zeros(size(img0));
vrImg_E(sub2ind(size(img0), round(vyESAC1_E), round(vxESAC1_E))) = 1;
vrImg_L(sub2ind(size(img0), round(vyESAC1_L), round(vxESAC1_L))) = 1;
H = fspecial('gaussian', [100 100], 20); %figure; imagesc(H);

vrImg_E = imfilter(vrImg_E, H); 
vrImg_L = imfilter(vrImg_L, H); 

figure; 
subplot 121; imshow(rgbmix(img0, vrImg_E)); title('Early');
subplot 122; imshow(rgbmix(img0, vrImg_L)); title('Late');
suptitle('ESAC>10 density');

%% ESAC amplitude as a distance from food