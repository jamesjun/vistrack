load D140330_Landmark;
load mlMask;

 %% plot trajectories
% 
% viAng = -1.15:.025:-1.0;
% figure;
% for i=1:numel(viAng)
%     subplot(6,1,i); 
%     imagesc(imadjust(imrotate(img0, viAng(i)))); 
%     title(num2str(viAng(i)));
%     axis([200 1400 610 640]);
%     colormap gray;
% end
% 
% %%
% mrImgF = fft2(double(img0));
% mrImgF = fftshift(real(mrImgF .* conj(mrImgF)));
% figure; imagesc(log(mrImgF(326:920, 504:1096)));

%%
angXaxis = -1.1590; %deg

vsPhase = {'E', 'L', 'P'};
mc = cell(3,2);
viAnimal = []; %plot all
nPlot = 4;
vcColor = 'bgrm';


figure;
% limInt = stretchlim(img0(mlMask));
for iPhase=1:3
    eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    
    
%     subplot(1,3,iPhase);  
    figure; 
    switch iPhase
        case 1
            img = img0;
            nPlot = 3;
            nTime = 2000; %10 sec
        case 2
            img = img0;
            nPlot = 3;
%             nTime = 100; %10 sec
        case 3
            img = img0_P;
            nPlot = 3;
            nTime = 2000; %10 sec
    end
    imshow(imadjust(img)); hold on;
    title(vsPhase{iPhase});
    viTrial = randperm(numel(vsTrialPool));   
    iPlot = 0;
    iTrialPlot = 1;
    while iPlot < nPlot
        S = vsTrialPool(viTrial(iTrialPlot));
        vrX = poolVecFromStruct(S, 'vrX', [], []);
        vrY = poolVecFromStruct(S, 'vrY', [], []); 
        if numel(vrX) >= nTime
            iPlot = iPlot + 1;
            if iPhase == 1 || iPhase == 2
                vrX = flipud(vrX(:));
                vrY = flipud(vrY(:));
            end
            if iPhase == 1 || iPhase == 3
                plot(vrX(1:nTime), vrY(1:nTime), vcColor(iPlot));             
            elseif iPhase == 2 %late. plot all
                plot(vrX, vrY, vcColor(iPlot));
            end
        end
        iTrialPlot = iTrialPlot + 1;
%         axis([504 1096 326 920]);
    end
end

%%
iTrialPlot = iTrialPlot + 1;
S = vsTrialPool_L(iTrialPlot);
vrX = poolVecFromStruct(S, 'vrX', [], []);
vrY = poolVecFromStruct(S, 'vrY', [], []); 
plot(vrX, vrY, 'b');


%% redo visit density, rotate images and square grid

%% Fig. 2A. Visit density map
pixpercm = 1053.28/(sqrt(2)*100);
rectCrop = [493 1083 312 902];
[RGB_E, mnVisit_E, mrTperV_E] = mapVisitCount(vsTrialPool_E);
[RGB_L, mnVisit_L, mrTperV_L] = mapVisitCount(vsTrialPool_L);
[RGB_P, mnVisit_P, mrTperV_P] = mapVisitCount(vsTrialPool_P);

figure; suptitle('Visit count'); 
subplot 131; imshow(RGB_E); title('Early'); axis(rectCrop);
subplot 132; imshow(RGB_L); title('Late');  axis(rectCrop);
subplot 133; imshow(RGB_P); title('Probe'); axis(rectCrop);


%% Fig. 5A. Back-swim locations
pixpercm = 1053.28/(sqrt(2)*100);
rectCrop = [493 1083 312 902];
[RGB_E, ~] = mapBackSwim(vsTrialPool_E);
[RGB_L, ~] = mapBackSwim(vsTrialPool_L);
[RGB_P, ~] = mapBackSwim(vsTrialPool_P);

figure; 
subplot 131; imshow(RGB_E); title('Early'); axis(rectCrop);
subplot 132; imshow(RGB_L); title('Late');  axis(rectCrop);
subplot 133; imshow(RGB_P); title('Probe'); axis(rectCrop);


%% Fig. 6A. Search time per grid
pixpercm = 1053.28/(sqrt(2)*100);
rectCrop = [493 1083 312 902];
[RGB_E, mrTperV_E] = mapSearchTime(vsTrialPool_E);
[RGB_L, mrTperV_L] = mapSearchTime(vsTrialPool_L);
[RGB_P, mrTperV_P] = mapSearchTime(vsTrialPool_P);

figure; suptitle('Visit count'); 
subplot 131; imshow(RGB_E); title('Early'); axis(rectCrop);
subplot 132; imshow(RGB_L); title('Late');  axis(rectCrop);
subplot 133; imshow(RGB_P); title('Probe'); axis(rectCrop);


%% %% crop and perform autocorr
figure;

nGrid = 20;
bootMean = @(x)[mean(x); bootci(1000, {@(y)mean(y), x})];
vcVisit = cell(3,1);

for iPhase = 1:3
% eval(sprintf('mr0 = mnVisit_%c;', vsPhase{iPhase})); 
eval(sprintf('mr0 = mrTperV_%c;', vsPhase{iPhase})); 
thresh = 1/exp(1);
img0ar = imadjust(imrotate(img0, -1.1590, 'nearest', 'crop'));
xyCr = rectCrop / nGrid;
mr = zeros(size(mr0));
viX = round(xyCr(1)):round(xyCr(2));
viY = round(xyCr(3)):round(xyCr(4));
mr(viY, viX) = mr0(viY, viX);
mr1 = imresize(mr, nGrid, 'nearest');

mrA = mr0(viY, viX); %active region cropped
cutoff = max(mrA(:)) * thresh;
vl = mrA >= cutoff;
ml1 = false(size(mr1)); ml1(mr1 >= cutoff) = 1;
vcVisit{iPhase} = double(vl(:));

subplot(3,2,1+(iPhase-1)*2);
mrVisit = uint8(mr1 / max(mr1(:)) * 255);
RGB = rgbmix(img0ar, mrVisit, mlMask);
imshow(RGB); axis(rectCrop);
subplot(3,2,2+(iPhase-1)*2);
imshow(maskRGB(RGB, ~ml1, .75)); axis(rectCrop);
title(vsPhase{iPhase})

end
figure; plotCellErrorbar(vcVisit, [], [], 1);
set(gca, 'XTickLabel', vsPhase);
ylabel('Fraction of search area');
axis([0 4 0 .3]); 
set(gca, 'YTick', 0:.05:.3);
title('Search area (>1/e) as a fraction of active zone');


% p value
[~, p] = kstest2(vcVisit{1}, [vcVisit{2}; vcVisit{3}])