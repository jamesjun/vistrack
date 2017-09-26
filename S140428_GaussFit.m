% S140428  2D gauss fitting

load D140330_Landmark;
load D140422_Controls;

nGrid = 20;
rectCrop = [493 1083 312 902]; %apply mask instead of square cut
xyFood0 = [789, 681];
xy0 = [787.0169  605.6858];
pixpercm = 1053.28/(sqrt(2)*100);

%% extract area to fit
vsTrialPool = vsTrialPool_P;
iAnimal = 1;

[RGB, mr0] = mapVisitCount(vsTrialPool, 0, iAnimal);

xyCr = round(rectCrop / nGrid);
mr1 = mr0(xyCr(3):xyCr(4), xyCr(1):xyCr(2));

figure; imagesc(mr1);

%
MdataSize = size(mr1,1);

[X,Y] = meshgrid(-MdataSize/2:MdataSize/2);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;

x0 = [1,MdataSize/2,MdataSize/4,MdataSize/2,MdataSize/4,0]; %Inital guess parameters

Z = D2GaussFunctionRot(x,xdata);
lb = [0,-MdataSize/2,0,-MdataSize/2,0,-pi/4];
ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2,pi/4];
[x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z,lb,ub);

%
%             amplitude of the 2D gaussian
%             centroid X of the 2D gaussian
%             centroid Y of the 2D gaussian
%             angle of rotation for the 2D gaussian
%             width X of the 2D gaussian
%             width Y of the 2D gaussian
[vmax, imax] = max(mr1(:));
vmin = min(mr1(:));
[iy, ix] = ind2sub(size(mr1), imax);
gof.StartPoint = [vmin, vmax, MdataSize/2, MdataSize/2, 0, 5, 5];
gof.Lower = [0, vmax/2, 1, 1, -180, 1, 1];
gof.Upper = [vmax/2, vmax*2, MdataSize, MdataSize, 180, 15, 20];
[fitresult, gof, Sfit] = Gauss2DRotFit(mr1, gof);
x0 = fitresult.c1;
y0 = fitresult.c2;
w0 = sqrt(fitresult.w1^2+fitresult.w2^2);
fprintf('x0=%f, y0=%f, w0=%f\n', x0,y0,w0); %translate back to the cm unit


% figure; imagesc(mr1); hold on; plot(x0,y0,'r.');

x0 = xyCr(1) * nGrid;
y0 = xyCr(3) * nGrid;
xyFit = ([fitresult.c1, fitresult.c2] + xyCr([1,3]) - [1.5 1.5]) * nGrid;

img0 = imadjust(vsTrialPool(1).img0);
figure; imshow(img0); 
hold on; plot(xyFood(1), xyFood(2), 'k.');
hold on; plot(xyFit(1), xyFit(2), 'r.');

xyFood = xyFood0 + vsTrialPool(1).xy0 - xy0;
DistErr_Mu = sqrt(sum((xyFit - xyFood).^2)) / pixpercm;
DistErr_SD = w0 * nGrid / pixpercm;

fprintf('iAnimal=%d, MU=%f cm, SD=%f cm\n', iAnimal, DistErr_Mu, DistErr_SD); %translate back to the cm unit

figure; imshow(img0); hold on; plot(xyFood(1), xyFood(2), 'k.');
%%
figure; imshow(RGB);
hold on; plot(xyFood(1), xyFood(2), 'k.');
hold on; plot(xyFit(1), xyFit(2), 'r.');
axis(rectCrop);
m = MdataSize;

mrX = reshape(Sfit.xData, m, m);
mrY = reshape(Sfit.yData, m, m);
mrZ = reshape(Sfit.zFit, m, m);
mrX = (mrX + xyCr(1) - 1.5) * nGrid;
mrY = (mrY + xyCr(3) - 1.5) * nGrid;
hold on; contour(mrX, mrY, mrZ);
axis square

figure; contour(mrX, mrY, mrZ); axis(rectCrop); axis square


%% Bar plot centre Mu


mrMu = [5.98 2.24 4.17 1.739; 8.338 6.448 10.038 4.826];
vrMu1 = mean(mrMu,2);
vrSd1 = std(mrMu,1,2)/2;

% shape plots, SD
mrSd = [29.32 24.03 30.298 22.199; 50.409 34.90 52.010 32.312];
vrMu2 = mean(mrSd,2);
vrSd2 = std(mrSd,1,2)/2;

figure;
subplot(1,2,1); hold on;
errorbar([1 2], vrMu1, vrSd1, 'k.');
plot(ones(1,4), mrMu(1,:), 'o');
plot(ones(1,4)*2, mrMu(2,:), 'o');
axis([.5 2.5 0 12]);
set(gca, 'XLim', [.5 2.5]);
set(gca, 'XTick', [1 2]); 
set(gca, 'XTickLabel', {'Stable LM', 'None LM'});
ylabel('Error Mean (cm)');
 
subplot(1,2,2); hold on;
errorbar([1 2], vrMu2, vrSd2, 'k.');
plot(ones(1,4), mrSd(1,:), 'o');
plot(ones(1,4)*2, mrSd(2,:), 'o');
axis([.5 2.5 0 60]);
set(gca, 'XLim', [.5 2.5]);
set(gca, 'XTick', [1 2]); 
set(gca, 'XTickLabel', {'Stable LM', 'None LM'});
ylabel('Error SD (cm)');
