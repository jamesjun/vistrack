% S140419 Heading direction analysis
load D140330_Landmark;



%% search-angle
warning off;
csCmd = {'abs(rad2deg(calcAngErr(RS.vrAh, RS.vrAf)))', '<|\theta_H - \theta_F|> (deg)', '@(x)mean(x)', '@(x)semcorr(x)'; ...
         'abs(rad2deg(calcAngErr(RS.vrAh, RS.vrAv)))', '<|\theta_H - \theta_V|> (deg)', '@(x)mean(x)', '@(x)semcorr(x)'; ...         
         'abs(rad2deg(RS.vrAt))', '<|\theta_T|> (deg)', '@(x)mean(x)', '@(x)semcorr(x)'};     
   
%      'rad2deg(calcAngErr(RS.vrAh, RS.vrAf))', 'SD \theta_H - \theta_F (deg)', '@(x)std(x)', '@(x)0'; ...         
         
iZone = 1; %plot all active
plotAll({vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}, csCmd, iZone);
suptitle('AZ');

%% angle correlation
warning off;
csCmd = {'IPI.vrCorrAhv_I', '<Corr(|\theta_H-\theta_V|, IPI)> (deg)', '@(x)mean(x)', '@(x)semcorr(x)'; ...
         'IPI.vrCorrAt_Ahv', '<Corr(|\theta_H-\theta_V|, \theta_T)> (deg)', '@(x)mean(x)', '@(x)semcorr(x)'; ...
         'IPI.vrCorrA_I', '<Corr(|\theta_T|, IPI)> (deg)', '@(x)mean(x)', '@(x)semcorr(x)'};
           
iZone = 1; %plot all active
plotAll({vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}, csCmd, iZone);
suptitle('AZ');


%%
angF = 129.8056; %from home to food
% S = poolTrials_RS(vsTrialPool_E);
% vrAv = normAng(rad2deg(S.vrAv) - angF);
% vrAe = normAng(rad2deg(S.vrAe));
% figure; ksdensity(abs(vrAv(S.vlZone)), 0:1:90, 'function', 'survivor', 'bandwidth', 1); 
% xlabel('|\theta_E1| (deg)'); set(gca, 'XLim', [0 90]);

% figure; hist((vrAe(S.vlZone), 0:1:90), -1:1:90, 'function', 'survivor', 'bandwidth', 1); 
% figure; hist(rad2deg(S.vrAv(S.vlZone)), -180:3:180);
% xlabel('|\theta_E2| (deg)'); set(gca, 'XLim', [-180 180]);
% figure; plot(S.vrAv, S.vrAf, '.');

figure;
strCmd = 'S.vrAt';
% strCmd = 'S.vrAh-deg2rad(angF)';
% strCmd = 'S.vrAh(S.vrDf<15)-S.vrAv(S.vrDf<15)';
% strCmd = 'calcAngErrSign(S.vrAv, S.vrAf)';
% strCmd = 'differentiate3(S.vrAv(:)-S.vrAh(:))*10';

% strSel = 'S.vrD1<3';
% strSel = 'S.vrD1<5 | S.vrD2<5 | S.vrD3<5 | S.vrD4<5';
strSel = 'S.vlZone';
% strSel = 'S.vrDf<15';

suptitle([strCmd, ' | ', strSel]);
csTrialPool = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P};
vcColor = 'rbg';
csPhase = {'E', 'L', 'P'};
for iPhase=1:3
    subplot(1,3,iPhase)
    S = poolTrials_RS(csTrialPool{iPhase});
    eval(sprintf('vr = %s;', strCmd));
    eval(sprintf('vl = %s;', strSel));
    rose(vr(vl), 360/5);
    h = get(gca, 'Children'); 
    set(h(1), 'Color', vcColor(iPhase));
    title(csPhase{iPhase});
end

%% swim direction. get tangent vector and plot histogram.

IPI = poolTrials_RS(vsTrialPool_L, 4);
vrTA = rad2deg(IPI.vrTA) - angF; %tangent angle

%if backward flip
% vrHA(IPI.vrV < 0) = vrHA(IPI.vrV < 0) + 180;
vrTA(vrTA < -180) = vrTA(vrTA < -180) + 360;
figure; ksdensity(vrTA);
set(get[
xlabel('Tagent angle');

%% error angle and IPIz

IPI = poolTrials_IPI(vsTrialPoooenl_L, 4);

vrTA = normAng(rad2deg(IPI.vrTA));

vrHA = rad2deg(IPI.vrHA); %head posture
vrHA(IPI.vrV < 0) = vrHA(IPI.vrV < 0) + 180;
vrHA = normAng(vrHA);
vrEA = normAng(vrHA - vrTA); %error angle

figure; plot(vrHA, vrTA, '.');
figure; ksdensity(vrEA); xlabel('Error Angle');

figure; plot(abs(vrEA), IPI.vrIz, '.'); %rho = -0.3114

figure; plot3(abs(vrEA), abs(IPI.vrV), IPI.vrIz, '.');

% kinematic model
calcVAF(abs(vrEA(:)), IPI.vrIz) %9.7%
calcVAF(abs(IPI.vrA), IPI.vrIz) %6% tail angle
calcVAF((IPI.vrV), IPI.vrIz) %18.0%
calcVAF([abs(vrEA(:)), abs(IPI.vrV(:)), abs(IPI.vrAc(:))], IPI.vrIz) %19.6% explained
calcVAF(abs(IPI.vrAc), IPI.vrIz) %5% explained
calcVAF([abs(vrEA(:)), abs(IPI.vrV(:)), abs(IPI.vrAc), abs(vrEA(:) .* IPI.vrV(:))], IPI.vrIz) %21.1%
calcVAF([abs(vrEA(:)), abs(IPI.vrV(:)), abs(IPI.vrAc), abs(vrEA(:) .* IPI.vrV(:)), abs(vrEA(:) .* IPI.vrAc(:))], IPI.vrIz) %21.4%


%% heading direction w.r.t the tangent of mid-point

warning off;
csCmd = {'abs(rad2deg(calcAngErr(RS.vrAh, RS.vrAf)))', '<|\theta_H - \theta_F|> (deg)', '@(x)mean(x)', '@(x)sem(x)'; ...
         'abs(rad2deg(calcAngErr(RS.vrAh, RS.vrAm)))', '<|\theta_H - \theta_M|> (deg)', '@(x)mean(x)', '@(x)sem(x)'; ...         
         'abs(rad2deg(RS.vrAt))', '<|\theta_T|> (deg)', '@(x)mean(x)', '@(x)sem(x)'};     
   
csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; csXlabel = {'Early', 'Late', 'Probe'};        
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csXlabel = {'Stable', 'None', 'Unstable'};
 
iZone = 2; %plot all active
plotAll(csTrials, csCmd, iZone, csXlabel);

%% angle correlation
warning off;

csCmd = {'RS.vrCorrV_At', 'V corr At', '@(x)mean(x)', '@(x)sem(x)'; ...
         'RS.vrCorrAe_At', 'Ae corr At', '@(x)mean(x)', '@(x)sem(x)'};

% csCmd = {'IPI.vrCorrIV', 'IPI corr V', '@(x)mean(x)', '@(x)semcorr(x)'};

% csCmd = {'IPI.vrCorrI1', 'Corr IPI vs \Delta\theta_H', '@(x)mean(x)', '@(x)semcorr(x)'; ...
%          'IPI.vrCorrI2', 'Corr IPI vs \Delta\theta_T', '@(x)mean(x)', '@(x)semcorr(x)'; ...
%          'IPI.vrCorrI3', 'Corr IPI vs |\theta_T|', '@(x)mean(x)', '@(x)semcorr(x)'; ...
%          'IPI.vrCorrAh_At', 'theta_oH corr theta_oT', '@(x)mean(x)', '@(x)semcorr(x)'};

% csCmd = {'IPI.vrCorrI_Vh', 'IPI corr Vh', '@(x)mean(x)', '@(x)semcorr(x)'};

csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; csXlabel = {'Early', 'Late', 'Probe'};        
% csTrials = {vsTrialPool_P, vsTrialPool_NP, vsTrialPool_WP}; csXlabel = {'Stable', 'None', 'Unstable'};

iZone = 1; %plot all active
plotAll(csTrials, csCmd, iZone, csXlabel);

%%
vx = -1:2:91;
iZone = 1;
RS = poolTrials_RS(csTrials{1},[],[],iZone);
vp1 = hist(rad2deg(RS.vrAe(RS.vlZ0)), vx);
RS = poolTrials_RS(csTrials{2},[],[],iZone);
vp2 = hist(rad2deg(RS.vrAe(RS.vlZ0)), vx);
RS = poolTrials_RS(csTrials{3},[],[],iZone);
vp3 = hist(rad2deg(RS.vrAe(RS.vlZ0)), vx);

figure; plot(vx, [vp1; vp2; vp3]);
legend({'E', 'L', 'P'});
% autocorr
% 
% vrAe = RS.vrAe;
% vrAe(vrAe>pi/2) = vrAe(vrAe>pi/2) -pi;
% figure; ksdensity(vrAe);

%% vector field

RS = poolTrials_RS(vsTrials{3});
X = RS.vrX;
Y = RS.vrY;


% vector alignment


%find grid average

figure; mapVelocityField(RS);


