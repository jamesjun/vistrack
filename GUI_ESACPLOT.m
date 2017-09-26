% GUI_ESACPLOT.m

%% Calculate 
LOADSETTINGS;

[EODR, TEOD, chName] = getSpike2Chan(handles.ADC, ADC_CH_EODR);
AMPL = getSpike2Chan(handles.ADC, ADC_CH_AMPL);
EODA = smoothFilter(differentiate5(EODR, .01), 5);
TLIM = handles.TC([1 end]);
IDXLIM = [];
IDXLIM(1) = find(TEOD > TLIM(1), 1, 'first');
IDXLIM(2) = find(TEOD < TLIM(2), 1, 'last');
IDXRNG = IDXLIM(1):IDXLIM(2);
TEOD = TEOD(IDXRNG);
EODR = EODR(IDXRNG);
EODA = EODA(IDXRNG);
[vrESAC, viESAC] = findESAC(EODA);
vtESAC = TEOD(viESAC);
vxESAC = interp1(handles.TC, filtPos(handles.XC(:,2), TRAJ_NFILT), TEOD(viESAC), 'spline');
vyESAC = interp1(handles.TC, filtPos(handles.YC(:,2), TRAJ_NFILT), TEOD(viESAC), 'spline');

% Kinematics calculation
[VELc, ACCc, ANGc, AVELc, XHc, YHc, TC] = calcVelocity(handles);
VEL = interp1(TC, VELc, TEOD, 'spline');
ACC = interp1(TC, ACCc, TEOD, 'spline');
ANG = interp1(TC, ANGc, TEOD, 'spline');
AVEL = interp1(TC, AVELc, TEOD, 'spline');
XH = interp1(TC, XHc, TEOD, 'spline');
YH = interp1(TC, YHc, TEOD, 'spline');
 
%% Focus on the top 5% of the ESAC
% triggered acceleration
pixpercm = 1053.28/(sqrt(2)*100);
vrRpix = sqrt((vxESAC - handles.xy0(1)).^2 + (vyESAC - handles.xy0(2)).^2);
vrRcm = vrRpix / pixpercm;
viWithin = vrRcm < 60;
viSelect = viESAC(viWithin);


%% Display
ccef = [];
vrLogESAC = log10(vrESAC(viWithin));

a=corrcoef(vrLogESAC, abs(VEL(viSelect))); ccef(1)=a(2);
a=corrcoef(vrLogESAC, abs(ACC(viSelect))); ccef(2)=a(2);
a=corrcoef(vrLogESAC, abs(ANG(viSelect))); ccef(3)=a(2);
a=corrcoef(vrLogESAC, abs(AVEL(viSelect))); ccef(4)=a(2);
a=corrcoef(abs(ANG(viSelect)), abs(ACC(viSelect))); ccef(5)=a(2);

% figure; 
% subplot 231; plot(vrLogESAC, abs(VEL(viSelect)), '.'); xlabel('Log10 ESAC'); ylabel('|VEL|'); 
%     title(sprintf('corr = %0.4f', ccef(1)));
% subplot 232; plot(vrLogESAC, abs(ACC(viSelect)), '.'); xlabel('Log10 ESAC'); ylabel('|ACC|'); 
%     title(sprintf('corr = %0.4f', ccef(2)));
% subplot 233; plot(vrLogESAC, abs(ANG(viSelect)), '.'); xlabel('Log10 ESAC'); ylabel('|ANG|'); 
%     title(sprintf('corr = %0.4f', ccef(3)));
% subplot 234; plot(vrLogESAC, abs(AVEL(viSelect)), '.'); xlabel('Log10 ESAC'); ylabel('|AVEL|'); 
%     title(sprintf('corr = %0.4f', ccef(4)));
% subplot 235; plot(abs(ANG(viSelect)), abs(ACC(viSelect)), '.'); xlabel('|ANG|'); ylabel('|ACC|'); 
%     title(sprintf('corr = %0.4f', ccef(5)));


% figure; 
% imshow(imadjust(handles.img0)); hold on;
% plot(vxESAC(viWithin), vyESAC(viWithin), 'r.');

%% 

figure;
subplot 221; plotBox(abs(ACC(viSelect)), vrESAC(viWithin)); ylabel('ESAC'); xlabel('|ACC|'); grid on;
subplot 222; plotBox(abs(ANG(viSelect)), vrESAC(viWithin)); ylabel('ESAC'); xlabel('|ANG|'); grid on;
subplot 223; plotBox(abs(AVEL(viSelect)), vrESAC(viWithin)); ylabel('ESAC'); xlabel('|AVEL|'); grid on;
subplot 224; plotBox(VEL(viESAC), vrESAC); ylabel('ESAC'); xlabel('VEL'); grid on;
% subplot 224; plotBox(sqrt(VEL(viSelect).^2+ANG(viSelect).^2), vrESAC(viWithin)); ylabel('ESAC'); xlabel('|ACC*ANG|'); grid on;

[~, dataID, ~] = fileparts(handles.vidFname);
suptitle(dataID);

%%
% collect top 5% esaccades within
vrESAC1 = vrESAC(viWithin);
viESAC1 = viESAC(viWithin);

viHigh = find(vrESAC1 >= quantile(vrESAC1, .95));
vrESAC2 = vrESAC1(viHigh);
viESAC2 = viESAC1(viHigh);

% triggered average of the acceleration
plotTrigAve((ACC), viESAC2, [-100 100]);  ylabel('ACC');
plotTrigAve((VEL), viESAC2, [-100 100]);  ylabel('VEL');

%% acceleration saccade triggered EODA
[vrASAC, viASAC] = findESAC(ACC);
% [vrASACn, viASACn] = findESAC(-ACC);

% select location within 60cm
% [vrASAC2, viASAC2] = selectWithinTop(XH, YH, vrASAC, viASAC, .95);

vxASAC = XH(viASAC);
vyASAC = YH(viASAC);
pixpercm = 1053.28/(sqrt(2)*100);
vrRpix = sqrt((vxASAC - handles.xy0(1)).^2 + (vyASAC - handles.xy0(2)).^2);
viWithin = find(vrRpix < (60 * pixpercm));
vrASAC1 = vrASAC(viWithin);
viASAC1 = viASAC(viWithin);
viHigh = find(vrASAC1 >= quantile(vrASAC1, .95));
vrASAC2 = vrASAC1(viHigh);
viASAC2 = viASAC1(viHigh);

plotTrigAve(abs(EODA), viASAC2, [-200 200]);  ylabel('EODA'); xlabel('Time from high ASAC');
% plotTrigAve((VEL), viASAC2, [-100 100]);  ylabel('VEL');
% figure; plot(ACC); hold on; plot(viASAC2, ACC(viASAC2), 'r.');

%%
% figure; 
% subplot 221; ksdensity(abs(ACC)); title('|ACC|');
% subplot 222; ksdensity(abs(ANG)); title('|ANG|');
% subplot 223; ksdensity(log(abs(ACC))); title('log |ACC|');
% subplot 224; ksdensity(log(abs(ANG))); title('log |ANG|');
