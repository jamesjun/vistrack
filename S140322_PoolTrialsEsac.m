% Pooled-trial analysis

dname = 'C:\expr\Landmark_Learning\';

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
    if iFish==1
        img0 = vsTrial{end}.img0; %save the background
    end
end

% S = makeStruct(TEOD, EODR, EODA, XH, YH, ...
%         VEL, ACC, ANG, AVEL, xy0, ...
%         viESAC, vrESAC, vtESAC, vxESAC, vyESAC, img0, dataID, pathLen_cm);

%% Plot Learning curve
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
    eval(sprintf('vrPathLen_%s = vrPathLen;', fishID));
    eval(sprintf('vrDuration_%s = vrDuration;', fishID));
end

figure;
subplot 221; bar(log10(vrPathLen_A),1); title('Fish A'); axis([0 50 1.5 4]); ylabel('log10 Dist');
subplot 222; bar(log10(vrPathLen_B),1); title('Fish B'); axis([0 50 1.5 4]); ylabel('log10 Dist');
subplot 223; bar(log10(vrPathLen_C),1); title('Fish C'); axis([0 50 1.5 4]); ylabel('log10 Dist');
subplot 224; bar(log10(vrPathLen_D),1); title('Fish D'); axis([0 50 1.5 4]); ylabel('log10 Dist');

figure;
subplot 221; bar(log10(vrDuration_A),1); title('Fish A'); axis([0 50 0 3]); ylabel('log10 Time');
subplot 222; bar(log10(vrDuration_B),1); title('Fish B'); axis([0 50 0 3]); ylabel('log10 Time');
subplot 223; bar(log10(vrDuration_C),1); title('Fish C'); axis([0 50 0 3]); ylabel('log10 Time');
subplot 224; bar(log10(vrDuration_D),1); title('Fish D'); axis([0 50 0 3]); ylabel('log10 Time');


%% pool 36 groups (4 sessions, 4 individuals = 16 data points)
vrPathMu = [];  vrPath_L = [];  vrPath_H = [];
vrDurMu = [];   vrDur_L = [];   vrDur_H = [];
mrPath = zeros(9, 16);
mrDur = zeros(9, 16);
for iSession=1:9
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
subplot 211; errorbar(1:9, vrPathMu, vrPath_L, vrPath_H); 
    ylabel('Path Length (cm)'); xlabel('Session #'); axis([.5 9.5 0 1500]);
subplot 212; errorbar(1:9, vrDurMu, vrDur_L, vrDur_H); 
    ylabel('Duration (s)'); xlabel('Session #'); axis([.5 9.5 0 150]);
suptitle('Landmark Learning Group, 4 animals pooled');


%% Select trials from early & late learning trials, keep mid 25-75%

mrPath_E = mrPath(1:3,:);
mrPath_L = mrPath(4:9,:);
pathQuant_E = quantile(mrPath_E(:), [.5 1]);
pathQuant_L = quantile(mrPath_L(:), [0 .5]);
mlPath_E = (mrPath_E >= pathQuant_E(1)) & (mrPath_E < pathQuant_E(2));
mlPath_L = (mrPath_L >= pathQuant_L(1)) & (mrPath_L < pathQuant_L(2));
mlPath = [mlPath_E; mlPath_L];

vrLogPath_E = log10(mrPath_E(:));
vrLogPath_L = log10(mrPath_L(:));
xi = 1.5:.1:4;
vnE = hist(vrLogPath_E, xi); 
vnL = hist(vrLogPath_L, xi);
vnE1 = hist(vrLogPath_E(mlPath_E), xi);
vnL1 = hist(vrLogPath_L(mlPath_L), xi);
figure; bar(xi, [vnE; vnL]', 1);
legend({'Early all', 'Late all'});
xlabel('log10 PathLen (cm)');
ylabel('# Obs.');

figure;
cdfplot(mrPath_E(:)); hold on;
cdfplot(mrPath_L(:)); hold on;

figure;
cdfplot(mrPath_E(mlPath_E)); hold on;
cdfplot(mrPath_L(mlPath_L)); hold on;


%% compare e-saccades distribution between early vs. late learning
appendArray = @(x,y)[x(:); y(:)];
%pool early
vrESAC_E = []; vrEsacEodr_E = []; vrEsacAcc_E = []; vrEsacAvel_E = []; vrEsacVel_E = [];
vrX_E = []; vrY_E = []; vrEODR_E = []; vrEODA_E = []; vrACC_E = []; vrVEL_E = []; vrANG_E = [];
vxESAC_E = []; vyESAC_E = [];
for iSession=1:3
    for iAnimal = 1:4
        eval(sprintf('vsTrial = vsTrial_%s;', vsFishID{iAnimal}));    
        for iTrial = 1:4
            if (mlPath(iSession, iTrial+(iAnimal-1)*4))
                S = vsTrial{iTrial+(iSession-1)*4};
                vlSelect = S.vlZone(S.viESAC);
                vrESAC_E = appendArray(vrESAC_E, S.vrESAC(vlSelect));
                vrEsacEodr_E = appendArray(vrEsacEodr_E, S.EODR(vlSelect));
                vrEsacAcc_E = appendArray(vrEsacAcc_E, S.ACC(vlSelect));
                vrEsacAvel_E = appendArray(vrEsacAvel_E, S.AVEL(vlSelect));
                vrEsacVel_E = appendArray(vrEsacVel_E, S.VEL(vlSelect));
                vxESAC_E = appendArray(vxESAC_E, S.vrX(vlSelect));
                vyESAC_E = appendArray(vyESAC_E, S.vrY(vlSelect));      
                
                vrX_E = appendArray(vrX_E, S.vrX(S.vlZone));
                vrY_E = appendArray(vrY_E, S.vrY(S.vlZone));
                vrEODR_E = appendArray(vrEODR_E, S.EODR(S.vlZone));
                vrEODA_E = appendArray(vrEODA_E, S.EODA(S.vlZone));
                vrACC_E = appendArray(vrACC_E, S.ACC(S.vlZone));
                vrVEL_E = appendArray(vrVEL_E, S.VEL(S.vlZone));
            end
        end
    end
end

%pool late
vrESAC_L = []; vrEsacEodr_L = []; vrEsacAcc_L = []; vrEsacAvel_L = []; vrEsacVel_L = [];
vrX_L = []; vrY_L = []; vrEODR_L = []; vrEODA_L = []; vrACC_L = []; vrVEL_L = []; vrANG_L = [];
vxESAC_L = []; vyESAC_L = [];
for iSession=4:9
    for iAnimal = 1:4
        eval(sprintf('vsTrial = vsTrial_%s;', vsFishID{iAnimal}));    
        for iTrial = 1:4
            if (mlPath(iSession, iTrial+(iAnimal-1)*4))
                S = vsTrial{iTrial+(iSession-1)*4};
                vlSelect = S.vlZone(S.viESAC); %location filtering
                vrESAC_L = appendArray(vrESAC_L, S.vrESAC(vlSelect));
                vrEsacEodr_L = appendArray(vrEsacEodr_L, S.EODR(vlSelect));
                vrEsacAcc_L = appendArray(vrEsacAcc_L, S.ACC(vlSelect));
                vrEsacAvel_L = appendArray(vrEsacAvel_L, S.AVEL(vlSelect));
                vrEsacVel_L = appendArray(vrEsacVel_L, S.VEL(vlSelect));
                vxESAC_L = appendArray(vxESAC_L, S.vrX(vlSelect));
                vyESAC_L = appendArray(vyESAC_L, S.vrY(vlSelect));
                
                vrX_L = appendArray(vrX_L, S.vrX(S.vlZone));
                vrY_L = appendArray(vrY_L, S.vrY(S.vlZone));
                vrEODR_L = appendArray(vrEODR_L, S.EODR(S.vlZone));
                vrEODA_L = appendArray(vrEODA_L, S.EODA(S.vlZone));
                vrACC_L = appendArray(vrACC_L, S.ACC(S.vlZone));
                vrVEL_L = appendArray(vrVEL_L, S.VEL(S.vlZone));               
            end
        end
    end
end


%% qq plot
figure; 
subplot 121; qqplot(vrESAC_E, vrESAC_L); axis([0 100 0 100]); axis square;
xlabel('E-Sac. Early (Hz/s)'); ylabel('E-Sac. Late (Hz/s)');
set(gca, {'XTick', 'YTick'}, {0:20:100, 0:20:100});

subplot 122; qqplot(vrESAC_E, vrESAC_L); axis([20 40 20 40]); axis square;
xlabel('E-Sac. Early (Hz/s)'); ylabel('E-Sac. Late (Hz/s)');
set(gca, {'XTick', 'YTick'}, {20:5:40, 20:5:40});

%%
threshESAC = 30;
calcProbEsac = @(y)[mean(y >= threshESAC); bootci(1000, {@(x)mean(x >= threshESAC), y})];
calcBootEsac = @(y)bootstrp(1000, @(x)mean(x >= threshESAC), y);
vec_E = calcProbEsac(vrESAC_E);
vec_L = calcProbEsac(vrESAC_L);

% vrESAC95_E = vrESAC_E(vrESAC_E >= quantile(vrESAC_E, .95));
% vrESAC95_L = vrESAC_L(vrESAC_L >= quantile(vrESAC_L, .95));
% calcMed = @(y)[median(y); bootci(1000, {@(x)median(x), y})];
% vec_E = calcMed(vrESAC_E);
% vec_L = calcMed(vrESAC_L);
figure; 
bar([1 2], [vec_E(1), vec_L(1)], .5);
hold on; 
errorbar([1 2], [vec_E(1), vec_L(1)], [vec_E(2), vec_L(2)], [vec_E(3), vec_L(3)], 'r.');
ylabel('Prob ESAC > 40');
set(gca, 'XTickLabel', {'Early', 'Late'})
[h p] = kstest2(calcBootEsac(vrESAC_E), calcBootEsac(vrESAC_L))

%Add errror bar
figure; ksdensity(vrESAC_E, 'function', 'survivor'); 
hold on; ksdensity(vrESAC_L, 'function', 'survivor'); 
legend({'Early', 'Late'});
xlabel('E-Saccades (Hz/s)');
ylabel('Survivor Function');
set(gca, 'YTick', 0:.01:.05);
axis([20 200 0 .05]);

figure; cdfplot(log10(vrESAC_L)); hold on; cdfplot(log10(vrESAC_E)); xlabel('log10 ESAC');
legend({'Late', 'Early'});

%% Plot the acceleration-triggered E-Sac

figure;
subplot 221; plotBox(vrESAC_E, abs(vrEsacAcc_E), [0 40 1000]); xlabel('ESAC'); ylabel('|Accel|'); title('Accel vs. E-Sac. Early');
subplot 222; plotBox(vrESAC_E, abs(vrEsacVel_E), [0 40 1000]); xlabel('ESAC'); ylabel('abs Vel.'); title('Accel vs. E-Sac. Early');
subplot 223; plotBox(vrESAC_L, abs(vrEsacAcc_L), [0 40 1000]); xlabel('ESAC'); ylabel('|Accel|'); title('Accel vs. E-Sac. Late');
subplot 224; plotBox(vrESAC_L, abs(vrEsacVel_L), [0 40 1000]); xlabel('ESAC'); ylabel('abs Vel.'); title('Accel vs. E-Sac. Late');

%% locations associated with ESAC: early vs. late. rotate and pool on the picture. 

figure;
subplot 121; plotBox(double(vrVEL_E>0), abs(vrEODA_E), [-1 .5 1]); 
ylabel('|EODA| (Hz/s)'); title('Early');  grid on; set(gca, 'XTick', [1 2]);
axis([.5 2.5 -20 20]);
set(gca, 'XTickLabel', {'Back', 'Forward'}); xlabel('Swim direction');

subplot 122; plotBox(double(vrVEL_L>0), abs(vrEODA_L), [-1 .5 1]); 
ylabel('|EODA| (Hz/s)'); title('Late'); grid on; set(gca, 'XTick', [1 2]);
axis([.5 2.5 -20 20]); 
set(gca, 'XTickLabel', {'Back', 'Forward'}); xlabel('Swim direction');

% [h p] = ttest2(vrEODR_L(vrVEL_L>0), vrEODR_L(vrVEL_L<0))

%% Triggered average of EODR before backswim

% backswim onset detection
calcProbBack = @(y)[mean(y < 0); bootci(1000, {@(x)mean(x < 0), y}, 'type', 'per')];
vec_E = calcProbBack(vrVEL_E);
vec_L = calcProbBack(vrVEL_L);

figure; 
bar([1 2], [vec_E(1), vec_L(1)], .5);
hold on; 
errorbar([1 2], [vec_E(1), vec_L(1)], [vec_E(2), vec_L(2)], [vec_E(3), vec_L(3)], 'r.');
ylabel('Prob. Backswim');
set(gca, 'XTickLabel', {'Early', 'Late'})
% [h p] = kstest2(calcBootEsac(vrESAC_E), calcBootEsac(vrESAC_L))


%% region-based analysis. location of e-saccades
threshESAC = 30;

vlESAC = vrESAC_E > threshESAC;
vxESAC1_E = vxESAC_E(vlESAC);
vyESAC1_E = vyESAC_E(vlESAC);

vlESAC = vrESAC_L > threshESAC;
vxESAC1_L = vxESAC_L(vlESAC);
vyESAC1_L = vyESAC_L(vlESAC);

figure; imshow(img0); hold on; title('Loc. E-Sacc > 40 during Early learning');
plot(vxESAC1_E, vyESAC1_E, 'b.');
plot(vxESAC1_L, vyESAC1_L, 'r.');