%% generate image mask
pixpercm = 1053.28/(sqrt(2)*100);
xy0 = [787 606];
xyFood = [789  681]; %food
xyLM1 = [966 418]; %landmark 1 (small square)
xyLM2 = [975 790]; %landmark 2 (large square)
xyLM3 = [604 799]; %landmark 3 (large circle)
xyLM4 = [600 428]; %landmark 4 (small circle)

%%
figure; imshow(imadjust(img0));
[h w] = size(img0);
[X Y] = meshgrid(1:w, 1:h);

mlMask = true(size(img0));

vrR0 = sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;
vrRfood = sqrt((X - xyFood(1)).^2 + (Y - xyFood(2)).^2) / pixpercm;

mlMask(vrR0 > 50 | vrRfood < 5) = 0;
mlMask(vrR0 > 60) = 0;

img0_m = img0; img0_m(~mlMask) = 0;
figure; imshow(imadjust(img0_m));

save mlMaskCentre50Food5 mlMask;

%%
calcDistCM = @(X, Y, xy0)sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;

figure; imshow(imadjust(img0));
[h w] = size(img0);
[X Y] = meshgrid(1:w, 1:h);

mlMask = false(size(img0));
mlMask(calcDistCM(X, Y, xyLM1) <= 10) = 1;
mlMask(calcDistCM(X, Y, xyLM2) <= 10) = 1;
mlMask(calcDistCM(X, Y, xyLM3) <= 10) = 1;
mlMask(calcDistCM(X, Y, xyLM4) <= 10) = 1;

img0_m = img0; img0_m(~mlMask) = 0;
figure; imshow(imadjust(img0_m));

save mlMaskLandmark10 mlMask;

%% LandmarkA
calcDistCM = @(X, Y, xy0)sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;

figure; imshow(imadjust(img0));
[h w] = size(img0);
[X Y] = meshgrid(1:w, 1:h);

mlMask = false(size(img0));
mlMask(calcDistCM(X, Y, xyLM4) <= 20) = 1;

img0_m = img0; img0_m(~mlMask) = 0;
figure; imshow(imadjust(img0_m));

save mlMaskLandmark20_D mlMask;

%% Food mask
calcDistCM = @(X, Y, xy0)sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;

figure; imshow(imadjust(img0));
[h w] = size(img0);
[X Y] = meshgrid(1:w, 1:h);

mlMask = false(size(img0));
mlMask(calcDistCM(X, Y, xyFood) <= 15) = 1;
mlMask(calcDistCM(X, Y, xyFood) <5) = 0;

img0_m = img0; img0_m(~mlMask) = 0;
figure; imshow(imadjust(img0_m));

save mlMaskFood5_15 mlMask;



%% LandmarkA
calcDistCM = @(X, Y, xy0)sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;

figure; imshow(imadjust(img0));
[h w] = size(img0);
[X Y] = meshgrid(1:w, 1:h);

mlMask = false(size(img0));
immask
% mlMask(calcDistCM(X, Y, xyLM1) <= 15) = 1;

img0_m = img0; img0_m(~mlMask) = 0;
figure; imshow(imadjust(img0_m));

save mlMaskLandmark_A mlMask;