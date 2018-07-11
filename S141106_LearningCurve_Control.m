% Pooled-trial analysis

%dname = 'C:\Users\James Jun\Desktop\New folder\';
% dname = '/Users/junj10/Desktop/None group learning/';
%dname = '/Users/junj10/Desktop/2013b_RandomNarrow_complete/';

mode = 'random';

switch lower(mode)
    case 'random'
%         dname = '/Users/junj10/Desktop/2013b_Random_complete/';
        dname = 'D:\Dropbox\Ms16_cue learning\vistrack\';
        outfname = 'D141026_RandGroup';
    case 'none'
        dname = '/Users/junj10/Desktop/None group learning/';
        outfname = 'D141026_NoneGroup';
    case 'stable'
        error('not implemented');
end
vsFnames = dir([dname '*_Track.mat']);
nFiles = numel(vsFnames);
vcFnames = {};
for i=1:nFiles
    vcFnames{i} = vsFnames(i).name;
end

% determine rotation and center
iData = 1;
S = load(fullfile(dname, vcFnames{iData}));
hFig = figure; 
imshow(imadjust(S.img0));
title('click (-50,0), (+50,0)cm');
set(gcf, 'Position', get(0, 'ScreenSize'));
[vrX, vrY] = ginput(2);
xy0 = [mean(vrX), mean(vrY)];
angXaxis = rad2deg(cart2pol(diff(vrX), diff(vrY))); %in rad
pixpercm = sqrt(diff(vrX)^2 + diff(vrY)^2) / 100;
hold on;
plot(vrX, vrY, 'r.');
plot(xy0(1), xy0(2), 'r.');
plot(S.xy0(1), S.xy0(2), 'go');
vcDisp = sprintf('iData: %d, ang: %0.3f deg, pixpercm: %0.3f, x0: %0.1f, y0: %0.1f', ...
    iData, angXaxis, pixpercm, xy0(1), xy0(2));
disp(vcDisp);
title(vcDisp);
uiwait(msgbox('press okay to continue'));
close(hFig); drawnow;

% Collect traces
% Rotate image?
warning off;
vsFishID = {'A', 'B', 'C', 'D'};
for iAnimal = 1:numel(vsFishID)
    fishID = vsFishID{iAnimal};
    vsTrial = [];
    csDataId = [];
    for iFile=1:nFiles
        fname = vcFnames{iFile};
        if fname(6) == '_' %exclude probe trials
            if strcmp(fname(4), fishID)
                vsTrial{end+1} = importTrial(fullfile(dname, fname), ...
                    'angXaxis', angXaxis, 'pixpercm', pixpercm);
                csDataId{end+1} = fname;
            end
        end
    end
    eval(sprintf('csDataId_%c = csDataId;', vsFishID{iAnimal}));
    eval(sprintf('vsTrial_%s = vsTrial;', fishID));
end

% S = makeStruct(TEOD, EODR, EODA, XH, YH, ...
%         VEL, ACC, ANG, AVEL, xy0, ...
%         viESAC, vrESAC, vtESAC, vxESAC, vyESAC, img0, dataID, pathLen_cm);
eval(sprintf('save %s vsTrial_*;', outfname));
fprintf('saved %s\n', outfname);


%%
mrColor = [0 0 1; 0 1 0; 1 0 1; 0 1 1; 1 0 0];
nSessions = 8;
nScores = 6; %store med, med_L, med_H, sd, sd_L, sd_H
vsFishID = {'A', 'B', 'C', 'D'};
csGroupFname = {'D140324_LandmarkGroup', 'D141026_NoneGroup', 'D141026_RandGroup'};
csGroup = {'StLm', 'NoLm', 'RnLm', 'RnLmCue', 'RnLm4'};
iGroupRand = 3;
csDataId = [];
nTrialsPerSession = 4;
nAnimals = numel(vsFishID);
nGroups = numel(csGroupFname);
trPathStat = zeros(nSessions, nGroups, nScores);
trDurStat = zeros(nSessions, nGroups, nScores);
trPath = zeros(nAnimals*nTrialsPerSession, nSessions, nGroups);
trDur = zeros(nAnimals*nTrialsPerSession, nSessions, nGroups);
cvPath = cell(nSessions, nGroups);
cvDur = cell(nSessions, nGroups);
for iGroup=1:numel(csGroupFname)
    vcPooledFname = csGroupFname{iGroup};
    load(vcPooledFname);

    for iAnimal = 1:numel(vsFishID)
        csDataId = [];
        animalID = vsFishID{iAnimal};
        eval(sprintf('vsTrial = vsTrial_%s;', animalID));
        vrPathLen = zeros(size(vsTrial));
        vrDuration = zeros(size(vsTrial));
        for iTrial = 1:numel(vsTrial)
            S = vsTrial{iTrial};
            try
            vrPathLen(iTrial) = S.pathLen_cm;
            vrDuration(iTrial) = diff(S.TEOD([1 end]));
            csDataId{end+1} = S.dataID;
            catch
                disp(S)
                vrPathLen(iTrial) = nan;
                vrDuration(iTrial) = nan;
            end
        end
        if iAnimal==1
            img0 = vsTrial{end}.img0; %save the background
        end   
        vrPathLen(isnan(vrPathLen)) = [];
        vrDuration(isnan(vrDuration)) = [];
        eval(sprintf('vrPathLen_%s = vrPathLen;', animalID));
        eval(sprintf('vrDuration_%s = vrDuration;', animalID));
        eval(sprintf('csDataId_%s_%s = csDataId;', csGroup{iGroup}, animalID));
    end

    for iSession=1:nSessions
        idxRng = (1:4) + (iSession-1)*4;
        vrPath = [vrPathLen_A(idxRng), vrPathLen_B(idxRng), vrPathLen_C(idxRng), vrPathLen_D(idxRng)];
        trPath(:, iSession, iGroup) = vrPath;
        trPathStat(iSession, iGroup, 1) = mean(vrPath);
        trPathStat(iSession, iGroup, 2) = median(vrPath);
        trPathStat(iSession, iGroup, 3) = std(vrPath);
        trPathStat(iSession, iGroup, 4) = iqr(vrPath);
        trPathStat(iSession, iGroup, 5) = quantile(vrPath, .25);
        trPathStat(iSession, iGroup, 6) = quantile(vrPath, .75);
        cvPath{iSession, iGroup} = vrPath(:);
        %trPathStat(iSession, iGroup, 5:6) = ...
        %    bootci(1000, {@(x)mean(x), vrPath});

        vrDur = [vrDuration_A(idxRng), vrDuration_B(idxRng), vrDuration_C(idxRng), vrDuration_D(idxRng)];
        trDur(:, iSession, iGroup) = vrDur;
        trDurStat(iSession, iGroup, 1) = mean(vrDur);
        trDurStat(iSession, iGroup, 2) = median(vrDur);
        trDurStat(iSession, iGroup, 3) = std(vrDur);
        trDurStat(iSession, iGroup, 4) = iqr(vrDur);
        trDurStat(iSession, iGroup, 5) = quantile(vrDur, .25);
        trDurStat(iSession, iGroup, 6) = quantile(vrDur, .75);
        cvDur{iSession, iGroup} = vrDur(:);
        %trDurStat(iSession, iGroup, 5:6) = ...
        %    bootci(1000, {@(x)mean(x), vrDur});
        
        if iGroup == iGroupRand
            for iGroup1=4:5
                switch iGroup1
                    case 4
                        vrPath1 = vrPath(1:end/2);
                        vrDur1 = vrDur(1:end/2);
                    case 5
                        vrPath1 = vrPath(end/2+1:end);
                        vrDur1 = vrDur(end/2+1:end);
                end
                trPathStat(iSession, iGroup1, 1) = mean(vrPath1);
                trPathStat(iSession, iGroup1, 2) = quantile(vrPath1, .5);
                trPathStat(iSession, iGroup1, 3) = std(vrPath1);
                trPathStat(iSession, iGroup1, 4) = iqr(vrPath1);
                trPathStat(iSession, iGroup1, 5) = quantile(vrPath1, .25);
                trPathStat(iSession, iGroup1, 6) = quantile(vrPath1, .75);
                cvPath{iSession, iGroup1} = vrPath1(:);
                %trPathStat(iSession, iGroup1, 5:6) = ...
                %    bootci(1000, {@(x)mean(x), vrPath1});

                trDurStat(iSession, iGroup1, 1) = mean(vrDur1);
                trDurStat(iSession, iGroup1, 2) = quantile(vrDur1, .5);
                trDurStat(iSession, iGroup1, 3) = std(vrDur1);
                trDurStat(iSession, iGroup1, 4) = iqr(vrDur1);
                trDurStat(iSession, iGroup1, 5) = quantile(vrDur1, .25);
                trDurStat(iSession, iGroup1, 6) = quantile(vrDur1, .75);
                cvDur{iSession, iGroup1} = vrDur1(:);
                %trDurStat(iSession, iGroup1, 5:6) = ...
                %    bootci(1000, {@(x)mean(x), vrDur1});
            end
        end
    end
end


%%
figure;
nRows = 2;
csPrefix = {'Mean', 'Median', 'SD', 'IQR'};
viPrefixPlot = [2 4];
viGroupPlot = [1 2 5];
nCols = numel(viPrefixPlot);
csGroup = {'LM', 'None', 'Rand', 'RandCue', 'Rand4'};
for iRow=1:nRows
    switch iRow
        case 1
            trData = trPathStat; 
            vcYlabel = 'Path Length (cm)'; 
            ytick = 0:1000:5000;
        case 2
            trData = trDurStat; 
            vcYlabel = 'Duration (s)'; 
            ytick = 0:100:500;
    end
    ylim = ytick([1 end]);
    
    for iCol=1:nCols
        iPrefix = viPrefixPlot(iCol);
        subplot(nRows, nCols, iCol + (iRow-1)*nCols); hold on;
        %for iGroup=1:size(trData,2)
        for iGroup = viGroupPlot
            errorbar(1:nSessions, trData(:,iGroup,iPrefix), 'color', mrColor(iGroup,:)); 
        end
        legend(csGroup{viGroupPlot});
        ylabel([csPrefix{iPrefix}, ' ', vcYlabel]); 
        xlabel('Session #');
        set(gca, 'XLim', [.5 (nSessions+.5)], 'XTick', 1:nSessions);
        set(gca, 'YLim', ylim, 'YTick', ytick);
    end
end

%% learning curve
figure;
csGroup = {'LM', 'None', 'Rand', 'RandCue', 'Rand4'};
viGroup = [1 2 5];
mrPath = trPathStat(:,viGroup,2);
mrPathL = trPathStat(:,viGroup,5);
mrPathH = trPathStat(:,viGroup,6);
mrDur = trDurStat(:,viGroup,2);
mrDurL = trDurStat(:,viGroup,5);
mrDurH = trDurStat(:,viGroup,6);
mrColor = [1 0 0; 0 1 0; 0 0 1];
for iRow=1:2
    switch iRow
        case 1
            mrY = mrPath / 100; 
            mrL = (mrY - mrPathL/100);
            mrH = (mrPathH/100 - mrY);
            vcYlabel = 'Path Length (m)'; 
            ytick = 0:10:60;
        case 2
            mrY = mrDur; 
            mrL = mrY - mrDurL;
            mrH = mrDurH - mrY;
            vcYlabel = 'Duration (s)'; 
            ytick = 0:100:600;
    end
    ylim = ytick([1 end]);
    
    subplot(nRows, 1, iRow);
    hold on;
    %bar(mrY, 1);
    for iGroup = 1:size(mrY,2)
        vrX = (1:nSessions) + iGroup/(size(mrY,2)+1) - .5;
        vrY = mrY(:,iGroup);
        h = bar(vrX, vrY, .2);
        set(h, 'FaceColor', mrColor(iGroup,:));
        set(h, 'EdgeColor', 'none');
        errorbar(vrX, vrY, mrL(:,iGroup), mrH(:,iGroup), 'color', mrColor(iGroup,:), 'LineStyle', 'none'); 
    end
    legend(csGroup{viGroup});
    ylabel(vcYlabel); 
    xlabel('Session #');
    set(gca, 'XLim', [.5 (nSessions+.5)], 'XTick', 1:nSessions);
    set(gca, 'YLim', ylim, 'YTick', ytick);
end


%% plot images
load('D141026_RandGroup');
trImg = [];
figure('Position', get(0, 'ScreenSize')); 
subplot(2,2,1);
for iTrial=1:numel(vsTrial)
    for iAnimal=1:nAnimals
        subplot(2,2,iAnimal);
        eval(sprintf('S = vsTrial_%c{iTrial};', vsFishID{iAnimal}));
        imshow(imadjust(S.img0));
        title(S.dataID);
    end
    pause(2);
end

%% compare across conditions
mrColor = [0 0 1; 0 1 0; 1 0 0];
cvTime = {1:2, 6:8};
csData = {'mean', 'std', 'mean3/4', 'std3/4', 'median', 'quantile1/4', 'quantile3/4'};
csType = {'Dur', 'Path'};
csTime = {'Early', 'Late'};
csGroup = {'StLm', 'NoLm', 'RnLm'};
trStat_Dur = zeros(numel(csData),numel(cvTime),5);
trStat_Path = zeros(numel(csData),numel(cvTime),5);
for iGroup=1:nGroups
    if iGroup == 3
        viSession = (size(trDur,1)/2 + 1):size(trDur,1);
    else
        viSession = 1:size(trDur,1);
    end
    for iTime = 1:2
        for iType = 1:2
            eval(sprintf('vr = tr%s(viSession,cvTime{iTime},iGroup);', csType{iType}));
            vr=sort(vr(:));
            n = round(numel(vr)/4);
            n1 = round(n/2); n2 = n - n1;
            vr1 = vr(n1+1:end-n2);
            vrStat = [mean(vr), std(vr), mean(vr1), std(vr1), median(vr), quantile(vr, 1/4), quantile(vr, 3/4)];
            eval(sprintf('trStat_%s(:,iTime,iGroup) = vrStat;', csType{iType}));
                 
        end% iType
    end %iTme
end %iGroup

%
csYlabel = {'Duration (s)', 'Path length (cm)'};
figure;
for iType=1:2
    subplot(2,1,iType); hold on;
    eval(sprintf('mrY = trStat_%s(1,:,:);\n', csType{iType}));
    eval(sprintf('mrE = trStat_%s(2,:,:);\n', csType{iType}));
    mrY = reshape(mrY(:), size(mrY,2), size(mrY,3));
    mrE = reshape(mrE(:), size(mrE,2), size(mrE,3));
    vcYlabel = csYlabel{iType};
    h = bar(mrY);
    for iGroup=1:nGroups
        vrXplot = (1:2) + .22 * (iGroup-2);
        set(h(iGroup), 'FaceColor', mrColor(iGroup,:), 'EdgeColor', 'none');
        errorbar(vrXplot, mrY(:,iGroup), zeros(size(mrE,1),1), ...
            mrE(:,iGroup), '.', 'Color', mrColor(iGroup,:));
    end
    legend(csGroup);
    set(gca, 'XTick', 1:2, 'XTickLabel', csTime);
    ylabel(vcYlabel);
end
title('mean & SD');

% pval
cvData = [];
for iType=1:2
    eval(sprintf('trData = tr%s;', csType{iType}));
    for iTime=1:2
        for iGroup=1:3
            if iGroup == 3
                viSession = (size(trDur,1)/2 + 1):size(trDur,1);
            else
                viSession = 1:size(trDur,1);
            end
            vr = trData(viSession,cvTime{iTime},iGroup);
            cvData{iGroup} = vr(:);
            [h, p] = lillietest(vr(:));
%             fprintf('%s:%s:%s, h: %f, p: %f\n', csType{iType}, csTime{iTime}, csGroup{iGroup}, h, p);
%             if h==0
%                 figure; hist(vr(:),10)
%             end
        end
        [h,pa]=ttest2(cvData{1}, cvData{2});
        [h,pb]=ttest2(cvData{1}, cvData{3});
        [h,pc]=ttest2(cvData{2}, cvData{3});
        fprintf('%s-%s\n', csType{iType}, csTime{iTime});
        fprintf('StLm-NoLm: %f\nStLm-RnLm: %f\nNoLm-RnLm: %f\n\n', pa, pb, pc);
    end
end

%% display trials selected
animalId = 'E';
nSessions = 8;
eval(sprintf('reshape(csDataId_RnLm_%c(1:4*%d), 4, %d)', animalId, nSessions, nSessions))

%% compare trials
%% Fig. 8B
viGroup = [1 2];
%Group: {'LM', 'None', 'Rand', 'RandCue', 'Rand4'};

cvTime = {1:2, 6:8};
quantLim = [1/8, 7/8];
mrColor = [1 0 0; 0 1 0; 0 0 1];
cvPathE = cell(size(viGroup));
cvPathL = cell(size(viGroup));

mrPathMu = zeros(3,2);
mrPathSd = zeros(3,2);

figure; hold on;
vcMode = 'Path'; %Path, Dur
for iMode = 1:2
    subplot(1,2,iMode);
    hold on;
    for iGroup=1:numel(viGroup)
        switch iMode 
            case 1
                cvData = cvPath;
                vcYlabel = 'Distance (m)';
                factor1 = 100;
            case 2
                cvData = cvDur;
                vcYlabel = 'Duration (s)';
                factor1 = 1;
        end

        cvPathE{iGroup} = quantFilt(cell2vec(cvData(cvTime{1}, viGroup(iGroup))), quantLim)/factor1;
        cvPathL{iGroup} = quantFilt(cell2vec(cvData(cvTime{2}, viGroup(iGroup))), quantLim)/factor1;
        mrPathMu(iGroup, 1) = nanmean(cvPathE{iGroup});
        mrPathMu(iGroup, 2) = nanmean(cvPathL{iGroup});
        mrPathSd(iGroup, 1) = bootSemMean(cvPathE{iGroup});
        mrPathSd(iGroup, 2) = bootSemMean(cvPathL{iGroup});
        vrX = (1:2) + iGroup/(numel(viGroup)+1) - .5;
        h = bar(vrX, mrPathMu(iGroup,:), .2);
        set(h, 'EdgeColor', 'None');
        set(h, 'FaceColor', mrColor(iGroup,:));
        errorbar(vrX, mrPathMu(iGroup,:), zeros(size(mrPathSd(iGroup,:))), mrPathSd(iGroup,:), 'k', 'LineStyle', 'none');
        [h,pa]=ttest2(cvPathE{iGroup}, cvPathL{iGroup});
        fprintf('Group %d: E vs L, p=%f\n', iGroup, pa);
    end
    kwtest_jjj([cvPathE cvPathL]);
   

    [h,p_E12]=ttest2(cvPathE{1}, cvPathE{2});
%     [h,p_E13]=ttest2(cvPathE{1}, cvPathE{3});
%     [h,p_E23]=ttest2(cvPathE{3}, cvPathE{2});
    [h,p_L12]=ttest2(cvPathL{1}, cvPathL{2});
%     [h,p_L13]=ttest2(cvPathL{1}, cvPathL{3});
%     [h,p_L23]=ttest2(cvPathL{3}, cvPathL{2});
    fprintf('Early: 1-2:%f, 1-3:%f, 2-3:%f\n', p_E12, p_E13, p_E23);
    fprintf('Late: 1-2:%f, 1-3:%f, 2-3:%f\n', p_L12, p_L13, p_L23);
    ylabel(vcYlabel);
    set(gca, 'XTick', 1:2); 
    set(gca, 'XTickLabel', {'Early', 'Late'});
end
% [Path]
% Group 1: E vs L, p=0.000007
% Group 2: E vs L, p=0.000033
% Group 3: E vs L, p=0.000001
% Early: 1-2:0.469637, 1-3:0.000007, 2-3:0.000039
% Late: 1-2:0.011141, 1-3:0.000142, 2-3:0.129105
% [Duration]
% Group 1: E vs L, p=0.000010
% Group 2: E vs L, p=0.000002
% Group 3: E vs L, p=0.000000
% Early: 1-2:0.705353, 1-3:0.000003, 2-3:0.000003
% Late: 1-2:0.070945, 1-3:0.000568, 2-3:0.067626