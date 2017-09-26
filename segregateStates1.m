function S = segregateStates1(AA, NRtime, thresh)
% S.TRISE, TFALL, TLENU, TLEND
% S.IDXU, IDXD (Index matched to NRtime
% assumes AA sampled at 100hz
% based on S89_moveMultichan

%% load data
% load('J:\2012.1W\Jan_prelesion\120120.mat');
% AA = V01_Ch701.values;
% NRtime = NRT0_120120;

%% compute

NBIN = 200; %3sec bin
SR = 100;

% thresh = -11; %determine from plotting the histogram
if size(AA,2) < 2
    VA = calcBinnedCV(AA, NBIN);
    VA(VA<0) = mean(VA);
    VAlog = log(VA)/log(10);

%     VA = calcBinnedSD(AA, NBIN);  
%     VAlog = log(VA)/log(10);
else
    VAlog = calcBinnedMovementMulti(AA, NBIN); %CV based, works better more specific    
%     VAlog = calcBinnedSDMulti(AA, NBIN);
%     VAlog = calcBinnedCVMultiMax(AA, NBIN); %CV based, works better more specific    
end

figure;
[PDF, XBIN, u] = ksdensity(VAlog, 'width', .1);
disp(sprintf('ks filter window: %f', u))
ksdensity(VAlog, 'width', .1);
% [N X] = ksdensity(VAlog);
% figure; plot(X, N);
% title(sprintf('Pulse Ampl CV Hist. NBIN = %d', NBIN));

%%
% thresh = -2.021;
if nargin < 3
    thresh = -1.5;
    [~, idxmin] = min(abs(XBIN - thresh));
    [idxmin] = findLocalMin(PDF, idxmin);
    thresh = XBIN(idxmin);
end
hold on; 
% threshLow = thresh * 1.1;
% threshHigh = thresh * 1;
threshLow = thresh - .1;
threshHigh = thresh + .1;
plot(thresh * [1 1], get(gca, 'YLim'), 'k');
plot(threshHigh * [1 1], get(gca, 'YLim'), 'r');
plot(threshLow * [1 1], get(gca, 'YLim'), 'b');

%%
% LA = VAlog > thresh;
% [LA, IRISE, IFALL, LENu, LENd] = cleanLogical(LA, 5, 6);


LA = threshHyst(VAlog, threshLow, threshHigh);
[LA, IRISE, IFALL, LENu, LENd] = parseLogical(LA, 1, 4);

IDXd = find(LA==0);
IDXu = find(LA==1);

figure; bar(IDXd, VAlog(IDXd), 1, 'b', 'EdgeColor',[0 0 1]);
hold on; bar(IDXu, VAlog(IDXu), 1, 'r', 'EdgeColor',[1 0 0]);
hold on; plot([1 numel(VAlog)], threshHigh * [1 1], 'r');
hold on; plot([1 numel(VAlog)], threshLow * [1 1], 'b');

% scale time
LA1 = repeatValues(LA, NBIN);

S.TRISE = (IRISE -.5) * NBIN / SR;
S.TFALL = (IFALL -.5)  * NBIN  / SR;
S.TLENU = LENu * NBIN  / SR;
S.TLEND = LENd * NBIN  / SR;
S.STATE = LA1;

%------------------------------------------------------------------------------
% separate distribution
if nargin > 1
    NRidx = round(NRtime * SR);
    ndata = numel(LA1);
    NRidx(NRidx > ndata) = [];
    S.IDXU = find(LA1(NRidx) == 1);
    S.IDXD = find(LA1(NRidx) == 0);
end

%% plot duration distribution set(gca, 'YScale', 'log');
figure;

subplot 221;
cdfplot(S.TLENU); set(gca, 'XScale', 'log'); 
title('Up state duration CDF');
xlabel('Time [sec]');
AX(1) = gca;

subplot 222;
cdfplot(S.TLEND); set(gca, 'XScale', 'log'); 
title('Down state duration CDF');
xlabel('Time [sec]');
AX(2) = gca;

subplot 223;
cdfplot(diff(S.TRISE)); set(gca, 'XScale', 'log'); 
title('rise-to-rise transition interval CDF');
xlabel('Time [sec]');
AX(3) = gca;

subplot 224;
cdfplot(diff(S.TRISE)); set(gca, 'XScale', 'log'); 
title('fall-to-fall transition interval CDF');
xlabel('Time [sec]');
AX(4) = gca;

linkaxes(AX, 'x');

%% show the state seperation with original Amplitude
% overlay state with the raw trace in bright yellow for up state and black for down state
if size(AA,2)>1
    AA = sum(AA,2);
end
T = (1:numel(AA))/SR;
figure; hold on; 
plot(T(1:4:end), AA(1:4:end), 'k');
nS = numel(LENu);
YLIM = get(gca, 'YLim');
height = abs(diff(YLIM));
for iS=1:nS    
    Tcoarse = (IRISE(iS) -.5) /SR*NBIN;
    POSVEC = [Tcoarse, YLIM(1), LENu(iS)/SR*NBIN, height];
    rectangle('Position',POSVEC,'LineWidth',1,'FaceColor',[1 0 0], 'EdgeColor', 'none');
end
plot(T(1:4:end), AA(1:4:end), 'k');
title('Amplitude state segregated');
axis tight;