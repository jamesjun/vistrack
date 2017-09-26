% create landmark masks

load D140330_Landmark;

%% create a mask
pixpercm = 1053.28/(sqrt(2)*100);

fig = figure; 
imshow(imadjust(img0))

%LM1. small square 2.216"
d1 = 2.216 * 2.54 * pixpercm; 
h1 = imrect(gca, [966-d1/2 418-d1/2 d1 d1]); %[x y w h]
mlMask1 = createMask(h1);

%LM2. Large square 3.545"
d2 = 3.545 * 2.54 * pixpercm; 
h2 = imrect(gca, [975-d2/2 790-d2/2 d2 d2]); %[x y w h]
mlMask2 = createMask(h2);

%LM3. large circle 4" dia
d3 = 4 * 2.54 * pixpercm; 
h3 = imellipse(gca, [604-d3/2 799-d3/2 d3 d3]); %[x y w h]
mlMask3 = createMask(h3);

%LM4. small circle 3" dia
d4 = 3 * 2.54 * pixpercm; 
h4 = imellipse(gca, [600-d4/2 428-d4/2 d4 d4]); %[x y w h]
mlMask4 = createMask(h4);

close(fig);

%% morphological operation

mlMask = mlMask1 | mlMask2 | mlMask3 | mlMask4;
% img1 = imadjust(img0); img1(~mlMask) = 0; 

mlMask0 = mlMask;
mlMask1 = imdilate(mlMask, strel('disk', round(1 * pixpercm))); %5 cm
mlMask2 = imdilate(mlMask, strel('disk', round(2 * pixpercm))); %5 cm
mlMask3 = imdilate(mlMask, strel('disk', round(3 * pixpercm))); %5 cm
mlMask4 = imdilate(mlMask, strel('disk', round(4 * pixpercm))); %5 cm
mlMask5 = imdilate(mlMask, strel('disk', round(5 * pixpercm))); %5 cm

img1 = imadjust(img0); img1(~mlMask5 | mlMask0) = 0; 
figure; imshow(img1);

%% build 5mm update
vsTrial = [vsTrialPool_E, vsTrialPool_L];
% xy0 = [789, 681]; %food
% xy0 = [966 418]; %LM1
% xy0 = [975 790]; %LM2
xy0 = [604 799]; %LM3
% xy0 = [600 428]; %LM4

vrEODA = poolVecFromStruct(vsTrial, 'EODA');
region = @(ml)getRegion(vsTrial, ml, xy0);

dr = .5;
vrR = 0:dr:8;
vrAsym = zeros(size(vrR));
for i=1:numel(vrAsym)    
    vrAsym(i) = calcDistAsym(vrEODA(region(getImageMask(img0, [0 dr] + vrR(i), 'LM3'))));
end
figure; bar(vrR+dr, vrAsym, 1); axis tight;
xlabel('Dist. from edge (cm)'); ylabel('EODAsym');
set(gca, 'XTick', 0:2:vrR(end)+dr);
title('LM3');

%% get IPI