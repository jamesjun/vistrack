% S140324 ESAC analysis



load D140324_LandmarkGroup;


%% Fig. 4D. EODA dist. during backswim
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