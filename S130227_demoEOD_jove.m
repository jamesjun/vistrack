% plot jve fig 2 E,F


fish = 'A';
lesion = 'pre'; %pre or post
dname = 'E:\SeagateExpr\2012.1W\fishA_prelesion_spont\';
dateid = '01-12_spont\02';

%% ------------------------------------------

%%
PARAM.thresh_AmplJump = .1;
PARAM.thresh_SD = -2;
% PARAM.thresh_DAmaxPost = .2;
PARAM.thresh_DAmaxPre= .05;
PARAM.thresh_Ampl0 = .1; %align amplitude
PARAM.TLIM = [-3 2];
PARAM.TLIMpre = [-9 0];
PARAM.TLIMpost = [0 1];
PARAM.nave = 4; %number of intervals to average with
PARAM.SR = 100;
PARAM.nwinAmpl = 10;
SR = PARAM.SR;

%%
trialType = [fish lesion 'Spont'];
dataid = ['12' dateid(1:2) dateid(4:5) '_' trialType];
load([dname dateid]);
[~, prefix] = fileparts(dateid);
eval(sprintf('PulseAmpl = V%s_Ch401.values;', prefix));
eval(sprintf('PulseTime = V%s_Ch401.times;', prefix));
% eval(sprintf('Rate = V%s_freeswim_000_Ch703.values;', prefix));
% eval(sprintf('DAmpl = V%s_freeswim_000_Ch704.values;', prefix));
% eval(sprintf('DRate = V%s_freeswim_000_Ch705.values;', dateid));
% eval(sprintf('Tmove = V%s_freeswim_000_Ch401.times;', prefix));

% load([dname dataid '_eod']);
% eval(sprintf('PulseTime = V%s_freeswim_000_Ch1.times;', prefix));
Ampl = abs(PulseAmpl);
Ampl = filtfilt(ones(1,PARAM.nwinAmpl), PARAM.nwinAmpl, Ampl);

%normalize
Ampl = Ampl / std(Ampl);
figure; plot(Ampl(100000:200000))

%% interpolate Rate
nPulsees = numel(PulseTime);
nave = PARAM.nave;
Rate1 = nave ./ (PulseTime(nave+1:end) - PulseTime(1:end-nave));
ILIM = round(nave/2) + [0 nPulsees - nave - 1];
Time1 = PulseTime(ILIM(1):ILIM(2));
% figure; plot(Time1(1:1000), Rate1(1:1000), '.-'); xlabel('Time (s)'); ylabel('Rate [Hz]');
Ampl1 = Ampl(ILIM(1):ILIM(2));

%% determine rate
TLIM = PulseTime([1 end]);
TLIM(2) = TLIM(1) + floor(diff(TLIM));

nt = diff(TLIM) * PARAM.SR;
Time = linspace(TLIM(1), TLIM(2), nt)';
Rate = interp1(Time1, Rate1, Time); %use linear interpolation
Ampl = interp1(Time1, Ampl1, Time); %use linear interpolation
DRate = (Rate(3:end) - Rate(1:end-2))/2 * (PARAM.SR);
DRate = [DRate(1); DRate(:); DRate(end)];
Time = Time - Time(1); %start from time zero. this is very important

%% segregate states
S = segregateStates1(Ampl1); %point process
Rup = Rate1(S.STATE);
Rdown = Rate1(~S.STATE);

%% plot IPI dist by states
figure; ksdensity(1./Rup, .01:.0001:.03);  set(gca, 'XLim', [.011 .021]);
figure; ksdensity(1./Rdown, .01:.0001:.03); set(gca, 'XLim', [.011 .021]); 

%% plot inst rate
figure; ksdensity(Rup, 30:.1:100);  set(gca, 'XLim', [45 85]);
figure; ksdensity(Rdown, 30:.1:100); set(gca, 'XLim', [45 85]); 

%% joint density

% plotJointDensity(X_int(~S.STATE), Y_int(~S.STATE), 200, [.01 .03 -.001 .001]); title('DI vs I. Down state. Pre');
% plotJointDensity(X_int(S.STATE), Y_int(S.STATE), 200, [.01 .03 -.001 .001]); title('DI vs I. Up state. Pre');
% plotJointDensity(X_int, Y_int, 200, [.01 .03 -.001 .001]); title('DI vs I. All states. Pre');

VA = calcBinnedCV(Ampl1, 200);
VA(VA<0) = mean(VA);
BinActivity = log(VA)/log(10);

BinIPI = calcBinnedMean(1./Rate1, 200);
BinRate = calcBinnedMean(Rate1, 200);

% N = plotJointDensity(BinActivity, BinRate, 200); title('Rate vs activity All states. Pre');


%% calc density costum
X = BinActivity;
Y = BinRate;
xbin = -3:.025:0;
ybin = 50:.5:85;
hfilt = [5 5];
nx = numel(xbin)-2;
ny = numel(ybin)-2;
N = zeros(ny, nx);
xplot = xbin(2:end-1);
yplot = ybin(2:end-1);
tic
for ix=1:nx
    IDX = find(X >= xbin(ix) & X < xbin(ix+2));
    for iy=1:ny
         N(iy,ix) = sum(Y(IDX) >= ybin(iy) & Y(IDX) < ybin(iy+2));
    end
end 
dx = mean(diff(xbin)) * 2;
dy = mean(diff(ybin)) * 2;
PDF = N / numel(X) / dx / dy; %prob norm
LogPdf = log(PDF)/log(10);
toc

% hGauss = fspecial('gaussian', hfilt, 2);
% N = imfilter(N, hGauss, 'replicate');
figure;
imagesc(LogPdf, 'XData', xplot, 'YData', yplot);
axis xy;
caxis([-3 -.5]);