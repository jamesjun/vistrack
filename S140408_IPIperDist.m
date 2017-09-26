% S140408_IPIperDist

load D140330_Landmark;
% create 4x6 figure tile for direct export to illustrator


%%

getIPIz = @(x)normIPI(differentiate3(poolVecFromStruct(x, 'vtEOD')));

getIPI = @(x)differentiate3(poolVecFromStruct(x, 'vtEOD'));
getDIPI = @(x)differentiate3(getIPI(x));

figure; plot(getIPI(vsTrialPool_E(1)), getDIPI(vsTrialPool_E(1)), '-');
getIPIn(vsTrialPool_E(1));

%% IPI normalization
figure; hold on;
for i=1:numel(vsTrialPool_E)
    vrIPIn = getDIPI(vsTrialPool_E(i));
    cdfplot(vrIPIn);
    h = get(gca, 'Children');
    set(h(1), 'Color', 'red');
end


for i=1:numel(vsTrialPool_L)
    vrIPIn = getDIPI(vsTrialPool_L(i));
    cdfplot(vrIPIn);
    h = get(gca, 'Children');
    set(h(1), 'Color', 'blue');
end

%% dL/L vs. dI. use LM3 circular

strRegion = 'CENTRE';
limDist = [0 50];
strLM = 'vrD4';
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

% dLdivL = @(x)differentiate3(x) ./ x;
dLdivL = @(x)x;
% dLdivL = @(x)x;

S_E = poolDIPIperDist(vsTrialPool_E, [], mlMask);
S_L = poolDIPIperDist(vsTrialPool_L, [], mlMask);
S_P = poolDIPIperDist(vsTrialPool_P, [], mlMask);

S = S_E;
vrDl = getfield(S, strLM);
vl1 = differentiate3(vrDl)<0 & vrDl < 5; % & S.vrV > 0; 
vl2 = differentiate3(vrDl)>0 & vrDl < 5; % & S.vrV > 0; 
figure; imshow(imadjust(img0)); hold on; 
plot(S.vrX(vl1), S.vrY(vl1), 'b.');
plot(S.vrX(vl2), S.vrY(vl2), 'r.');
disp([mean(S.vrD(vl1)), mean(S.vrD(vl2))])

%% plot vrD as a function of distance
S = poolDIPIperDist(vsTrialPool_L, [], mlMask);
% strTitle = 'Distance to the edge of landmarks';
strTitle = 'Distance to the centre of food';
% strTitle = 'Distance to the edge of food';

fun2 = @(x)1./  x;

fun1 = @(I, D, lim)(I(D >= lim(1) & D < lim(2)));

dx = .5;
vrX = (0:dx/2:30);
vrY = zeros(size(vrX));
vrL = zeros(size(vrX));
vrH = zeros(size(vrX));
for i=1:numel(vrX)
    dlim = vrX(i) + [-1 1] * dx/2;

%     vr = [fun1(S.vrD, S.vrD1, dlim); fun1(S.vrD, S.vrD2, dlim); ...
%           fun1(S.vrD, S.vrD3, dlim); fun1(S.vrD, S.vrD4, dlim)]; %LM
       
%     vr = [fun1(S.vrD, S.vrD1, dlim); fun1(S.vrD, S.vrD2, dlim); ...
%           fun1(S.vrD, S.vrD3, dlim); fun1(S.vrD, S.vrD4, dlim); ...
%           fun1(S.vrD, S.vrDf, dlim)]; %all

%     vr = [fun1(S.vrD, S.vrDf, dlim)]; %food only
    
%     vr = [fun1(S.vrD, S.vrDf, dlim)]; %food only, probe

    vr = [fun1(S.vrD, S.vrDf+1, dlim)]; %food only, probe
    
    [vrY(i), vrL(i), vrH(i)] = bootSEM(vr, fun2);
end
hold on; plotError(vrX, vrY, vrL, vrH, 'r');
set(gca, 'XLim', vrX([1 end]));
title(strTitle);
% xlabel('Distance (cm)');
% ylabel('IPI/cm');
% set(gca, 'YLim', [0 2]);  set(gca, 'YTick', 0:.5:2);
% set(gca, 'XTick', vrX(1):2:vrX(end));


%% approach vs. depart from the landmarks, comparison of DIPI

strTitle = 'Dist to the edge of landmarks';

fun1 = @(x)1./  x / 10;
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
mc = cell(3,2);
% vrCorrD = [];
for i=1:3
    eval(sprintf('S = poolDIPIperDist(vsTrialPool_%s, [], mlMask);', vsPhase{i}));
    mc{i,1} = S.vrD(vlLM_app(S));
    mc{i,2} = S.vrD(vlLM_dep(S));
%     vrCorrD = [vrCorrD(:); S.vrCorrD(:)];
    [~,p] = kstest2(mc{i,1}, mc{i,2});
    disp(p);
end

figure; plotCellErrorbar(mc, '', fun1, 1);
title('IPI/mm during landmark approach (blue) or depart (red)');
xlabel('Trial Phase');
ylabel('IPI/mm');
set(gca, 'YTick', 0:.25:1);
set(gca, {'XTick', 'XTickLabel'}, {1:3, vsPhase});
axis([.5 3.5 0 1]);


%% IPI/mm Forward vs. backswim
strRegion = 'CENTRE';
limDist = [0 50];
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

strTitle = ['Forward vs. backward swim central region R<50cm'];

fun1 = @(x)1./  x / 10;
dlim = [0 3];

vlLM = @(S)(S.vrD1 >= dlim(1) & S.vrD1 < dlim(2)) |...
           (S.vrD2 >= dlim(1) & S.vrD2 < dlim(2)) |...
           (S.vrD3 >= dlim(1) & S.vrD3 < dlim(2)) |...
           (S.vrD4 >= dlim(1) & S.vrD4 < dlim(2));
       
vlFb = @(S)(S.vrDf >= 3 & S.vrDf < 14);

vlFa = @(S)(S.vrDf >= dlim(1) & S.vrDf < dlim(2));   

vlAll = @(S)true(size(S.vrD));

vsPhase = {'E', 'L', 'P'};
mc = cell(3,2);
for i=1:3
    eval(sprintf('S = poolDIPIperDist(vsTrialPool_%s, [], mlMask);', vsPhase{i}));
    mc{i,1} = S.vrD(S.vrV > 0 & vlAll(S));
    mc{i,2} = S.vrD(S.vrV < 0 & vlAll(S));
end

figure; plotCellErrorbar(mc, '', fun1, 1);
% axis([0 4 0 1.5]);
ylabel('IPI/mm');
% set(gca, {'YLim', 'YTick'}, {[0 1.25], 0:.25:1.25});

xlabel('Region');
set(gca, {'XTick', 'XTickLabel'}, {1:3, {'Early', 'Late', 'Probe'}});
% set(gca, {'XTick', 'XTickLabel'}, {1:4, {'R<50', 'LM<3', 'F3~10', 'F<3'}});

title(strTitle);


%% IPI/mm by region
strRegion = 'CENTRE';
limDist = [0 50];
[mlMask, regionStr] = getImageMask(img0, limDist, strRegion);

strTitle = ['IPI/mm by regions & learning-associated changes'];

fun1 = @(x)1./  x / 10;
dlim = [0 3];

vlLM = @(S)(S.vrD1 >= dlim(1) & S.vrD1 < dlim(2)) |...
           (S.vrD2 >= dlim(1) & S.vrD2 < dlim(2)) |...
           (S.vrD3 >= dlim(1) & S.vrD3 < dlim(2)) |...
           (S.vrD4 >= dlim(1) & S.vrD4 < dlim(2));
       
vlFb = @(S)(S.vrDf >= 3 & S.vrDf < 14);

vlFa = @(S)(S.vrDf >= dlim(1) & S.vrDf < dlim(2));   

vlAll = @(S)true(size(S.vrD));

vsPhase = {'E', 'L', 'P'};
mc = cell(4,3);
for i=1:3
    eval(sprintf('S = poolDIPIperDist(vsTrialPool_%s, [], mlMask);', vsPhase{i}));
    mc{1,i} = S.vrD(vlAll(S));
    mc{2,i} = S.vrD(vlLM(S));
    mc{3,i} = S.vrD(vlFb(S));
    mc{4,i} = S.vrD(vlFa(S));
end

figure; plotCellErrorbar(mc, '', fun1, 1);
% axis([0 4 0 1.5]);
ylabel('IPI/mm');
set(gca, {'YLim', 'YTick'}, {[0 1], 0:.25:1});

xlabel('Region');
set(gca, {'XTick', 'XTickLabel'}, {1:4, {'R<50', 'LM<3', 'F3~14', 'F<3'}});

title(strTitle);


%% grid bitmap image


%% EODAs
% [mnVisit1, mnVisit] = calcGridStats(vsTrialPool_E, img0, 'EODAs');
titleStr = 'IPI/mm, Late trials';
vsPhase = {'E', 'L', 'P'};
rectCrop = [493 1083 312 902];

figure; 
for iPhase = 1:3
    subplot(1,3,iPhase);
    eval(sprintf('vsTrial = vsTrialPool_%s;', vsPhase{iPhase}));
    S = poolDIPIperDist(vsTrial, []);
    [mnVisit, mnVisit1, mnI] = mapIPIperMM(S);
    imshow(mnI);
    axis(rectCrop);
    title(vsPhase(iPhase));
end

%% xcorr2

% mlMask1 = imresize(mlMask, 1/gridSize);
% ylim = round([351,860]/gridSize); xlim = round([522, 1053]/gridSize);
ylim = round([459,760]/gridSize); xlim = round([636,943]/gridSize);
mnVisitC = mnVisit(ylim(1):ylim(2), xlim(1):xlim(2));
figure; 
mrXcorr = xcorr2(mnVisit);
imagesc(mrXcorr); %imCropMask(mnVisit, mlMask1)));
axis square; set(gca, {'XTick', 'YTick'}, {[],[]});
caxis auto;
title(titleStr);
[a b] = calcXcorr2Max(mrXcorr, .6);
disp([a b]);