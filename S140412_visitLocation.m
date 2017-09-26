% visit location plot

load D140330_Landmark;

%% sum landmark areas
d = 3;
r1 = 2.216*2.54/2; %*1.1222; %cm, radius
r2 = 3.545*2.54/2; %*1.1222; %cm, radius
r3 = 4*2.54/2; %cm, radius
r4 = 3*2.54/2; %cm, radius
A_LM = pi*((r1+d)^2 - r1^2) + pi*((r2+d)^2 - r2^2) + (pi*d^2 + d*4*2*r3) + (pi*d^2 + d*4*2*r4);
A_LMa = pi*(r1^2) + pi*(r2^2) + (2*r3)^2 + (2*r4)^2;
A_Z = 6400 - A_LMa;
A_F = pi*4^2;
A_Fa = pi*15^2;

%% all animals pooled
vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,1);
for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_location(vsTrialPool, []);
%     vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
    vl = S.vrDf <= 14;
    cvLM{iPhase} = double(vl(:)) / sum(S.vlZone) * numel(vl); % to be pooled    
end


figure;
plotCellErrorbar(cvLM);
set(gca, {'XTick', 'XTickLabel'}, {1:3, vsPhase});
set(gca, 'XLim', [.5 3.5]);
set(gca, {'YLim', 'YTick'}, {[0 .5], 0:.1:1});
% [~, p] = kstest2_jjj(cvLM{1}, cvLM{2})
% [~, p] = kstest2_jjj(cvLM{2}, cvLM{3})
title('Fraction of time near food (R<4cm)');
hold on; plot(get(gca, 'XLim'), A_Fa/A_Z*[1 1], 'k');

%% analyze per animal

vsPhase = {'E', 'L', 'P'};
cvLM = cell(4,3);
for iAnimal = 1:4
    for iPhase = 1:3
        eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
        S = poolTrials_location(vsTrialPool, iAnimal);
        vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
%         vl = S.vrDf <= 3;
        cvLM{iAnimal,iPhase} = double(vl(:)) / sum(S.vlZone) * numel(vl); % to be pooled    
    end
end
figure;
plotCellErrorbar(cvLM);
title('Fraction of time near landmarks (<3cm)');
set(gca, {'XTick', 'XTickLabel'}, {1:4, {'A', 'B', 'C', 'D'}});
xlabel('');

%% show region definitions
rectCrop = [493 1083 312 902];
img0a = imadjust(img0);
[mlMask1, regionStr] = getImageMask(img0, [0 3], 'LM*F');
[mlMask2, regionStr] = getImageMask(img0, [0 0], 'LM*F');
mlMask = mlMask1 & ~mlMask2;
img0a = imrotate(img0a, -1.1590, 'nearest', 'crop');
% img0a(~mlMask) = img0a(~mlMask) * .5;
figure; imshow(img0a);
axis(rectCrop);



%% back swim prob. all animals pooled
vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,1);
for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_location(vsTrialPool, []);
    
    % Region
    vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
%     vl = S.vrDf < 3;
%     vl = logical(S.vlZone);
    
    cvLM{iPhase} = double(S.vrV(vl) < 0);
end


figure;
plotCellErrorbar(cvLM);
set(gca, {'XTick', 'XTickLabel'}, {1:3, vsPhase});
set(gca, 'XLim', [.5 3.5]);
set(gca, 'YLim', [0 .16]);
set(gca, 'YTick', 0:.04:.16);
% [~, p] = kstest2_jjj(cvLM{1}, cvLM{2})
% [~, p] = kstest2_jjj(cvLM{2}, cvLM{3})
title('Back-swim prob in landmarks<3');
% hold on; plot(get(gca, 'XLim'), A_Fa/A_Z*[1 1], 'k');

%% analyze per animal

vsPhase = {'E', 'L', 'P'};
cvLM = cell(4,3);
for iAnimal = 1:4
    for iPhase = 1:3
        eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
%         S = poolTrials_location(vsTrialPool, iAnimal);
        S = poolTrials_IPI(vsTrialPool, iAnimal);
%         vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
        vl = S.vrDf <= 14;
%         vl = logical(S.vlZone);
        
%         cvLM{iAnimal, iPhase} = double(S.vrV(vl) < 0);
        cvLM{iAnimal, iPhase} = abs(S.vrDA(vl));
    end
end
figure;
plotCellErrorbar(cvLM);
title('Back-swim prob food <15cm');
set(gca, {'XTick', 'XTickLabel'}, {1:4, {'A', 'B', 'C', 'D'}});
xlabel('');


%% tail angle all animals pooled
vsPhase = {'E', 'L', 'P'};
vcColor = 'rbg';
cvLM = cell(3,1);
figure; hold on;
for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_location(vsTrialPool, []);
%     S = poolTrials_IPI(vsTrialPool, []);
    
    % Region
    vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
%     vl = S.vrDf < 14;
%     vl = logical(S.vlZone);
% %     vl = S.vrV > 0;
    vrZ = abs(S.vrA(vl)) * 180 / pi;
%     ksdensity(vrZ, 'function', 'survivor'); 
%     h=get(gca, 'Children'); 
%     set(h(1), 'color', vcColor(iPhase));
    
    vrAl = differentiate3(differentiate3(S.vrA));
    vrW = (-vrAl(vl) ./ S.vrA(vl));
    vrW = sqrt(vrW(vrW>0));
%     disp(mean(vrW));    
    ksdensity(vrW, 0:.005:.5, 'function', 'survivor'); 
    h=get(gca, 'Children'); 
    set(h(1), 'color', vcColor(iPhase));
    
    cvLM{iPhase} = vrW;
end
legend({'E', 'L', 'P'});
xlabel('|Tail-Ang| [deg]'); 
axis([0 100 0 1]);

figure;
plotCellErrorbar(cvLM);
set(gca, {'XTick', 'XTickLabel'}, {1:3, vsPhase});
set(gca, 'XLim', [.5 3.5]);
% set(gca, 'YLim', [0 .16]);
% set(gca, 'YTick', 0:.04:.16);
% [~, p] = kstest2_jjj(cvLM{1}, cvLM{2})
% [~, p] = kstest2_jjj(cvLM{2}, cvLM{3})
title('Fc<15');
ylabel('|Tail-Ang| [deg]');
% hold on; plot(get(gca, 'XLim'), A_Fa/A_Z*[1 1], 'k');

%% tail bending angle
vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,1);
figure;
for iPhase = 1:3
    subplot(1,3,iPhase);
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    [RGB, mnVisit] = mapTailBending(vsTrialPool);
    imshow(RGB);
    title(vsPhase{iPhase});
    axis(rectCrop);
end

%% change of Tail-bending per IPI

vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,1);
figure;
for iPhase = 1:3
    subplot(1,3,iPhase); hold on;
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_IPI(vsTrialPool, []);
    for iLM=1:4
        eval(sprintf('vrDi = S.vrD%d;', iLM));
        vrV = differentiate3(vrDi);
%         vrV = vrV(vrDi<3);
%         vrV = log(abs(vrV(vrDi < 3) ./ vrDi(vrDi < 3)));
        plot(abs(vrV), differentiate3(S.vrD(vrDi<3)), '.');
    end
end
ylabel('Diff D/ipi'); xlabel('|V|');


%% tail angle all animals pooled
vsPhase = {'E', 'L', 'P'};
cvLM = cell(4,3);
figure;
for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_IPI(vsTrialPool, []);
    
    % Region
%     
%     vl = S.vrDf < 3;
      
%     vl = logical(S.vlZone);
%     vl = S.vrV > 0;
%     subplot(1,3,iPhase);
%     plot(S.vrA, S.vrDA, '-');
%     cvLM{iPhase} = abs(S.vrDA(vl));
%     title(vsPhase{iPhase});
%     vl = vl & abs(S.vrDA) < .005;
    
    for iRegion = 1:4
        switch iRegion
            case 1
                vl = S.vlZone;
            case 2
                vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
            case 3
                vl = S.vrDf < 14;
            case 4
                vl = S.vrDf < 3;
        end        
        cvLM{iRegion, iPhase} = abs(S.vrDA(vl)) * 180 / pi;
%         cvLM{iRegion, iPhase} = S.vrD(vl)*10;
    end
end


figure;
plotCellErrorbar(cvLM, '', @(x)1./x,1);
set(gca, 'XTickLabel', {'AZ', 'LM<3', 'Fc<15', 'F<3'});
ylabel('IPI/deg');
set(gca, 'YLim', [0 1.25]); set(gca, 'YTick', 0:.25:1.25);
% set(gca, 'YLim', [0 2]); set(gca, 'YTick', 0:.5:2);
% [~, p] = kstest2_jjj(cvLM{1}, cvLM{2})
% [~, p] = kstest2_jjj(cvLM{2}, cvLM{3})
title('IPI/deg, tail-bending angle');
% title('IPI/mm, distance between EOD');
% hold on; plot(get(gca, 'XLim'), A_Fa/A_Z*[1 1], 'k');


%% tail bending angle
vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,1);
figure;
for iPhase = 1:3
    subplot(1,3,iPhase);
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_IPI(vsTrialPool, []);
    [RGB, mnVisit] = mapTailBending(S);
    imshow(RGB);
    title(vsPhase{iPhase});
    axis(rectCrop);
end


%% IPI/angle forward vs. backward
vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,2);

for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_IPI(vsTrialPool, []);

%     vl = S.vlZone;
    vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone   
    
    cvLM{iPhase, 1} = abs(S.vrDA(vl & S.vrV > 0))*180/pi;
    cvLM{iPhase, 2} = abs(S.vrDA(vl & S.vrV < 0))*180/pi;
end


figure;
plotCellErrorbar(cvLM, '', @(x)1./x, 1);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
ylabel('IPI/deg');
set(gca, 'YLim', [0 1.5]); set(gca, 'YTick', 0:.5:2);
% set(gca, 'YLim', [0 2]); set(gca, 'YTick', 0:.5:2);
% [~, p] = kstest2_jjj(cvLM{1}, cvLM{2})
% [~, p] = kstest2_jjj(cvLM{2}, cvLM{3})
% title('IPI/deg, tail-bending angle');
title('IPI/deg, Forward vs. Backward, F<3');
% hold on; plot(get(gca, 'XLim'), A_Fa/A_Z*[1 1], 'k');

%% approach vs depart
dlim = [0 3];
vlLM_dep = @(S)(S.vrD1 >= dlim(1) & S.vrD1 < dlim(2) & differentiate3(S.vrD1) > 0) |...
               (S.vrD2 >= dlim(1) & S.vrD2 < dlim(2) & differentiate3(S.vrD2) > 0) |...
               (S.vrD3 >= dlim(1) & S.vrD3 < dlim(2) & differentiate3(S.vrD3) > 0) |...
               (S.vrD4 >= dlim(1) & S.vrD4 < dlim(2) & differentiate3(S.vrD4) > 0);
         
vlLM_app = @(S)(S.vrD1 >= dlim(1) & S.vrD1 < dlim(2) & differentiate3(S.vrD1) < 0) |...
               (S.vrD2 >= dlim(1) & S.vrD2 < dlim(2) & differentiate3(S.vrD2) < 0) |...
               (S.vrD3 >= dlim(1) & S.vrD3 < dlim(2) & differentiate3(S.vrD3) < 0) |...
               (S.vrD4 >= dlim(1) & S.vrD4 < dlim(2) & differentiate3(S.vrD4) < 0);

vsPhase = {'E', 'L', 'P'};
cvLM = cell(3,2);

for iPhase = 1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolTrials_IPI(vsTrialPool, []);

%     vl = S.vlZone;
%     vl = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone   
    vl_app = vlLM_app(S);
    vl_dep = vlLM_dep(S);
    cvLM{iPhase, 1} = abs(S.vrDA(vl_app))*180/pi;
    cvLM{iPhase, 2} = abs(S.vrDA(vl_dep))*180/pi;
end


figure;
plotCellErrorbar(cvLM, '', @(x)1./x, 1);
set(gca, 'XTickLabel', {'Early', 'Late', 'Probe'});
ylabel('IPI/deg');
set(gca, 'YLim', [0 1.5]); set(gca, 'YTick', 0:.5:2);
% set(gca, 'YLim', [0 2]); set(gca, 'YTick', 0:.5:2);
% [~, p] = kstest2_jjj(cvLM{1}, cvLM{2})
% [~, p] = kstest2_jjj(cvLM{2}, cvLM{3})
% title('IPI/deg, tail-bending angle');
title('IPI/deg, Approach vs. Depart, LM<3');
% hold on; plot(get(gca, 'XLim'), A_Fa/A_Z*[1 1], 'k');


%% DIPI stats

plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, 'vrI', @(x)calcCorrTau(x));
plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, 'vrDI', @(x)calcCorrTau(x));
plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, 'vrD', @(x)calcCorrTau(x));
plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, 'vrDI', @(x)std(x));

%%
pixpercm = 7.4478;

% MOV stats
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'abs(RS.vrV) / 7.4478', [], '<|Vel|> (cm/s)'); set(AX(1), 'YLim', [0 15]);
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'double(RS.vrV<0)', [], 'Prob. Backswim'); set(AX(1), 'YLim', [0 .5]);
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'RS.vrV / 7.4478', @(x)std(x), 'SD Vel (deg/s)');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    '(RS.vrV / 7.4478).^2', @(x)mean(x), '<Vel^2> (deg/s)^2');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    '(rad2deg(RS.vrA))', @(x)std(x), 'SD T.Ang (deg)');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'rad2deg(abs(RS.vrAV))', @(x)mean(x), '|T.Avel| (deg/s)');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'rad2deg(RS.vrAV)', @(x)std(x), 'SD T.Avel (deg/s)');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'rad2deg(abs(RS.vrHAV))', @(x)mean(x), '<|H.Avel|> (deg/s)');
[AX, AX1, cvZ] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'rad2deg(RS.vrHA)', @(x)std(x), 'SD H.Ang (deg)');

% EOD stats
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrI*1000', @(x)mean(x), '<IPI> (ms)'); set(AX(1), 'YLim', [12.5 14.5]);
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrI*1000', @(x)skewness(x), 'Skewness IPI (ms)');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrDI*1000', @(x)std(x), 'SD D.IPI (ms)'); %set(gca, 'YLim', [12.5 14.5]);
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrDI*1000', @(x)skewness(x), 'Skewness D.IPI (ms)'); %set(gca, 'YLim', [12.5 14.5]);
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    '(IPI.vrDI*1000).^2', @(x)mean(x), '<D.IPI^2> (ms^2)');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrDI * 1000', @(x)mean(x), '<DI> (ms)');

% Both
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrD', @(x)mean(x*10), '<mm/IPI>');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrDA', @(x)mean(abs(rad2deg(x))), '<T.Ang deg/IPI>'); 
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrDHA', @(x)mean(abs(rad2deg(x))), '<H.Ang deg/IPI>'); 
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrCorrIV(:)', @(x)mean(x), 'Corr IPI-Vel');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrCorrIA(:)', @(x)mean(x), 'Corr IPI-A.Ang');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrCorrIHA(:)', @(x)mean(x), 'Corr IPI-H.Ang');

[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrCorrID(:)', @(x)mean(x), 'Corr mm/IPI-IPI');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrCorrVD(:)', @(x)mean(x), 'Corr mm/IPI-Vel');
[AX, AX1] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrCorrIV(:)', @(x)mean(x), 'IPI-Vel');

%% Variance accounted for
sem = @(x)std(x)/numel(x);
% vsTrialPool = [vsTrialPool_E, vsTrialPool_L, vsTrialPool_P];
% vsTrialPool = [vsTrialPool_L];
vsTrialPool = [vsTrialPool_E];

IPI_A = poolTrials_IPI(vsTrialPool, 1);
IPI_B = poolTrials_IPI(vsTrialPool, 2);
IPI_C = poolTrials_IPI(vsTrialPool, 3);
IPI_D = poolTrials_IPI(vsTrialPool, 4);

mrMean = [mean(IPI_A.vrCorrID), mean(IPI_A.vrCorrVD), mean(IPI_A.vrCorrIVD); ...
          mean(IPI_B.vrCorrID), mean(IPI_B.vrCorrVD), mean(IPI_B.vrCorrIVD); ...
          mean(IPI_C.vrCorrID), mean(IPI_C.vrCorrVD), mean(IPI_C.vrCorrIVD); ...
          mean(IPI_D.vrCorrID), mean(IPI_D.vrCorrVD), mean(IPI_D.vrCorrIVD)];
mrSem = [sem(IPI_A.vrCorrID), sem(IPI_A.vrCorrVD), sem(IPI_A.vrCorrIVD); ...
          sem(IPI_B.vrCorrID), sem(IPI_B.vrCorrVD), sem(IPI_B.vrCorrIVD); ...
          sem(IPI_C.vrCorrID), sem(IPI_C.vrCorrVD), sem(IPI_C.vrCorrIVD); ...
          sem(IPI_D.vrCorrID), sem(IPI_D.vrCorrVD), sem(IPI_D.vrCorrIVD)]; 
figure; plotBarError(mrMean, mrSem, {'m', 'c', 'k'}, {'A', 'B', 'C', 'D'});
ylabel('Corr. mm/IPI vs. {IPI, Speed, Combined}');

mrMean = [mean(IPI_A.vrVafID), mean(IPI_A.vrVafVD), mean(IPI_A.vrVafIVD); ...
          mean(IPI_B.vrVafID), mean(IPI_B.vrVafVD), mean(IPI_B.vrVafIVD); ...
          mean(IPI_C.vrVafID), mean(IPI_C.vrVafVD), mean(IPI_C.vrVafIVD); ...
          mean(IPI_D.vrVafID), mean(IPI_D.vrVafVD), mean(IPI_D.vrVafIVD)];
mrSem = [sem(IPI_A.vrVafID), sem(IPI_A.vrVafVD), sem(IPI_A.vrVafIVD); ...
         sem(IPI_B.vrVafID), sem(IPI_B.vrVafVD), sem(IPI_B.vrVafIVD); ...
         sem(IPI_C.vrVafID), sem(IPI_C.vrVafVD), sem(IPI_C.vrVafIVD); ...
         sem(IPI_D.vrVafID), sem(IPI_D.vrVafVD), sem(IPI_D.vrVafIVD)]; 
figure; plotBarError(mrMean, mrSem, {'m', 'c', 'k'}, {'A', 'B', 'C', 'D'});
ylabel('VAF. mm/IPI vs. {IPI, Speed, Combined}');

%%
[AX, AX1, cvZ] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrVafIV(:)', [], 'VAF Vel by IPI'); axis([.5 1.5 0 1])
[AX, AX1, cvZ] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, ...
    'IPI.vrVafVI(:)', [], 'VAF IPI by Vel'); axis([.5 1.5 0 1])


%% cdf plot swim speed, IPI, and IPI/mm. check for normalization scheme
% option = 3; %vrI, vrD
pixpercm = 7.4478;
% if ~exist('cvZmu', 'var')
    cvZmu = cell(4,3);
    cvZsd = cell(4,3);
% end
% vsTrial = [vsTrialPool_L, vsTrialPool_L, vsTrialPool_P];
vsTrial = [vsTrialPool_L]; %stable phase only
vcColor = 'krmgb';
for option=1:3
    figure; hold on;
    for iAnimal=1:4
        for iTrial=1:numel(vsTrial)
            if vsTrial(iTrial).iAnimal == iAnimal
                S = poolTrials_IPI(vsTrial(iTrial));
%                 S1 = poolTrials_RS(vsTrial(iTrial));
                switch option
                    case 1, vrZ = S.vrI * 1000;
                    case 2, vrZ = S.vrD * 10;
                    case 3, vrZ = abs(S.vrV / pixpercm);
                end
                vrZ = vrZ(S.vlZone);
                cvZmu{iAnimal, option} = [cvZmu{iAnimal, option}; mean(vrZ)];
                cvZsd{iAnimal, option} = [cvZsd{iAnimal, option}; std(vrZ)];
                ksdensity(vrZ, 'Function', 'survivor');
                h = get(gca, 'Children');
                set(h(1), 'Color', vcColor(iAnimal+1));
            end
        end
    end
    switch option
        case 1, xlabel('IPI (ms)'); axis([10 16 0 1]);
        case 2, xlabel('Dist./IPI (mm)'); axis([0 5 0 1]);
        case 3, xlabel('Speed (cm/s)'); %axis([0 5 0 1]);
    end
end

%%
figure; hold on;
vrX = [];   vrY = [];
ix = 2;
iy = 3;
vrCrr = [];
for iAnimal = 1:4
    vrX1 = cvZmu{iAnimal,ix};
    vrY1 = cvZmu{iAnimal,iy};
    c = corrcoef(vrX1, vrY1);
    vrCrr(end+1) = c(2);
    vrX = [vrX; vrX1(:)];
    vrY = [vrY; vrY1(:)];
    plot(vrX, vrY, '.');
%     , '.', 'color', vcColor(iAnimal));
end
c = corrcoef(vrX, vrY);
vrCrr(end+1) = c(2);
title(sprintf('Mu corr = %f', vrCrr(end)));
switch ix
    case 1, xlabel('IPI (ms)');
    case 2, xlabel('Dist/IPI (mm)');
    case 3, xlabel('Speed (cm/s)');
end
switch iy
    case 1, ylabel('IPI (ms)');
    case 2, ylabel('Dist/IPI (mm)');
    case 3, ylabel('Speed (cm/s)');
end
disp(vrCrr)

%% animal specific distribution
iy = 3;
vrX = [];
vrY = [];
for iAnimal = 1:4
    vrY1 = cvZmu{iAnimal,iy};    
    vrY = [vrY; vrY1(:)];
    vrX = [vrX; iAnimal * ones(size(vrY1))];
end

csAnimal = {'A', 'B', 'C', 'D'};
x_ordinal = ordinal(vrX, csAnimal, [], .5:1:5);

figure; 
boxplot(vrY, x_ordinal); %plot except the rest, superimpose

switch iy
    case 1, ylabel('IPI (ms)');
    case 2, ylabel('Dist/IPI (mm)');
    case 3, ylabel('Speed (cm/s)');
end