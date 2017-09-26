function S = poolTrials_IPI(vsTrial, iAnimal, lim, iZone)
% iAnimal=0 for all
%Distance to landmark and IPI

pixpercm = 1053.28/(sqrt(2)*100);
%landmark locations
xyf = [789, 681];    rf = 1; %cm, radius
xy0 = vsTrial(1).xy0;
xy1 = [966, 418];    r1 = 2.216*2.54/2; %*1.1222; %cm, radius
xy2 = [975, 790];    r2 = 3.545*2.54/2; %*1.1222; %cm, radius
xy3 = [604, 799];    r3 = 4*2.54/2; %cm, radius
xy4 = [600, 428];    r4 = 3*2.54/2; %cm, radius

if nargin < 2, iAnimal = []; end
if nargin < 3, lim = []; end
if nargin < 4, iZone = 1; end

%Select animal-specific trials
if ~isempty(iAnimal) && iAnimal > 0
    viAnimal = poolVecFromStruct(vsTrial, 'iAnimal');
    vsTrial = vsTrial(viAnimal == iAnimal);
end

calcD0 = @(x,y)sqrt((x-xy0(1)).^2 + (y-xy0(2)).^2) / pixpercm;        
calcD1 = @(x,y)dist2square((x-xy1(1))/pixpercm, (y-xy1(2))/pixpercm, r1);
calcD2 = @(x,y)dist2square((x-xy2(1))/pixpercm, (y-xy2(2))/pixpercm, r2);
calcD3 = @(x,y)sqrt((x-xy3(1)).^2 + (y-xy3(2)).^2) / pixpercm - r3;        
calcD4 = @(x,y)sqrt((x-xy4(1)).^2 + (y-xy4(2)).^2) / pixpercm - r4;        
calcDf = @(x,y)sqrt((x-xyf(1)).^2 + (y-xyf(2)).^2) / pixpercm - rf;        
calcSpeed = @(x,y)sqrt(differentiate5(x).^2 + differentiate5(y).^2) * 100 / pixpercm;
calcAng = @(x,y)cart2pol(differentiate5(x), differentiate5(y)); %tangent angle in rad
interpAng = @(t,y,t1)interp1(t, unwrap(y), t1, 'spline');

S = [];
S.img0 = vsTrial(1).img0;
S.xy0 = vsTrial(1).xy0; 

for iTrial=1:numel(vsTrial)
    Strial = vsTrial(iTrial);
    
    %trim
    Strial = trimField(Strial, lim, ...
        'TEOD', 'vrX', 'vrY', 'TEOD', 'VEL', 'ANG', 'HANG', 'AVEL', 'HAVEL');    
    vrTr = Strial.TEOD;
    vtEOD = Strial.vtEOD;
    vtEOD(vtEOD < vrTr(1) | vtEOD > vrTr(end)) = [];    
    
    % interpolate the location of the EOD pulse
    vrXr = Strial.vrX;
    vrYr = Strial.vrY;
    vrLr = cumsum(sqrt(differentiate3(vrXr).^2 + differentiate3(vrYr).^2));
    vrX = interp1(vrTr, vrXr, vtEOD, 'spline');
    vrY = interp1(vrTr, vrYr, vtEOD, 'spline');
    vrL = interp1(vrTr, vrLr, vtEOD, 'spline')  / pixpercm; %pathlen
    vrV = interp1(vrTr, Strial.VEL, vtEOD, 'spline') / pixpercm;
    vrAcc = interp1(vrTr, Strial.ACC, vtEOD, 'spline') / pixpercm;        
    vrVh = interp1(vrTr, calcSpeed(vrXr, vrYr), vtEOD, 'spline') / pixpercm;
    
    %angle
    mrXr = Strial.mrX;
    mrYr = Strial.mrY;
    vrAvr = unwrap(cart2pol(differentiate3(Strial.vrX), differentiate3(Strial.vrY))); %tangential angle
    vrWvr = differentiate5(vrAvr);
    vrAv = interp1(vrTr, vrAvr, vtEOD, 'spline');    
    vrAh = interpAng(vrTr, Strial.HANG, vtEOD);
    vrAt = interp1(vrTr, Strial.ANG, vtEOD, 'spline');    
    vrWv = interp1(vrTr, vrWvr, vtEOD, 'spline');
    vrWh = interp1(vrTr, Strial.HAVEL, vtEOD, 'spline');
    vrWt = interp1(vrTr, Strial.AVEL, vtEOD, 'spline');
    vrALt = interp1(vrTr, differentiate5(Strial.AVEL), vtEOD, 'spline'); %ang accel
    vrAf = cart2pol(xyf(1)-vrX, xyf(2)-vrY);
    vrAhf = calcAngErr(vrAh, vrAf);
    vrAhv = calcAngErr(vrAh, vrAv);
    vrAvf = calcAngErr(vrAv, vrAf);
    vrAm = interpAng(vrTr, calcAng(mrXr(:,4), mrYr(:,4)), vtEOD);
    vrAdt = interpAng(vrTr, cart2pol(mrXr(:,5)-mrXr(:,6), mrYr(:,5)-mrYr(:,6)), vtEOD);
    vrAdh = interpAng(vrTr, cart2pol(mrXr(:,2)-mrXr(:,3), mrYr(:,2)-mrYr(:,3)), vtEOD);
    vrAdt_Am = calcAngErr(vrAdt, vrAm);
    vrAdh_Am = calcAngErr(vrAdh, vrAm);
    
    vrAhm = calcAngErr(vrAh, vrAm);
    
    %per IPI
%     vrD = sqrt(differentiate3(vrX).^2 + differentiate3(vrY).^2) / pixpercm;      
    vrD = differentiate3(vrL);
    vrDAv = differentiate3(vrAv);
    vrDAh = differentiate3(vrAh);
    vrDAt = differentiate3(vrAt);    
    vrAv = wrap(vrAv);
    vrAh = wrap(vrAh);
%     vrDA = interp1(vrTr, poolVecFromStruct(Strial, 'AVEL'), vtEOD, 'spline');
    vrI = differentiate3(vtEOD);
    vrIz = zscore(vrI);
    vrIq = qscore(vrI);
    vrDI = differentiate3(vrI);
    vrD0 = calcD0(vrX, vrY);
    vrD1 = calcD1(vrX, vrY);
    vrD2 = calcD2(vrX, vrY);
    vrD3 = calcD3(vrX, vrY);
    vrD4 = calcD4(vrX, vrY);
    vrDf = calcDf(vrX, vrY);
    vlZone = isZone(vrX, vrY, Strial.xy0);
    vrPmm = 1/mean(vrD(vlZone)*10);
    vrPmm_F15 = 1/mean(vrD(vrDf<15)*10);
    
    %DIPI E-sac in zone
    viTs = findDIsac(vrDI);    
    vrDsac = diff(vrL(viTs));
    vrIsac = diff(vtEOD(viTs));
    vrXsac = vrX(viTs(2:end));
    vrYsac = vrY(viTs(2:end));
    viTs0 = viTs(2:end);
    vlZone_s = isZone(vrX(viTs0), vrY(viTs0), Strial.xy0);
%     vrDsac = vrDsac0(vlZone_s);
%     vrIsac = vrIsac0(vlZone_s);
    try
    vlZone_s = vrD1(viTs0)<3 | vrD2(viTs0)<3 | vrD3(viTs0)<3 | vrD4(viTs0)<3;
    vrDsac_LM = vrDsac(vlZone_s);
    vrIsac_LM = vrIsac(vlZone_s);
    vlZone_s = vrDf(viTs0)<15;
    vrDsac_F15 = vrDsac(vlZone_s);
    vrIsac_F15 = vrIsac(vlZone_s);    
    catch
        disp('err');
    end
    
%     vrCorrV = calcCorrTau(vrV);
%     vrCorrD = calcCorrTau(vrD);
%     vrCorrI = calcCorrTau(vrI);
    
    S1 = makeStruct(vrD0, vrD1, vrD2, vrD3, vrD4, vrDf, vlZone);
    vlZ0 = getZone(S1, iZone);
    vrCorrI1 = calcCorr(vrI, vrAdh_Am, vlZ0);
    vrCorrI2 = calcCorr(vrI, vrAdt_Am, vlZ0);
    vrCorrI3 = calcCorr(vrI, abs(vrAt), vlZ0);
    vrCorrIV = calcCorr(vrI, abs(vrV), vlZ0);
    vrCorrI_Vh = calcCorr(vrI, abs(vrVh), vlZ0);
    
    vrCorrAh_At = calcCorr(abs(vrAdh), abs(vrAdt), vlZ0);

    S = appendField(S, vrX, vrY, vrV, vrAcc, ...
            vrI, vrDI, vrIz, vrIq, ...
            vrD, vrDAt, vrDAh, vrDAv, vrPmm, vrPmm_F15, ...
            vrDsac, vrIsac, vrDsac_LM, vrIsac_LM, vrDsac_F15, vrIsac_F15, ...
            vrD0, vrD1, vrD2, vrD3, vrD4, vrDf, vlZone, ...
            vrAh, vrAt, vrAv, vrWh, vrWt, vrWv, vrALt, ...
            vrAf, vrAhm, vrAhv, ...
            vrCorrI1, vrCorrI2, vrCorrI3, vrCorrAh_At, ...
            vrCorrIV, vrCorrI_Vh, vrXsac, vrYsac);            
    
%     vrCorrAdt_Am = calcCorr(vrAdt, vrAm, vlZ0);
%     vrCorrAdh_Am = calcCorr(vrAdh, vrAm, vlZ0);
%         
%     vrVafA_VHA = calcVAF(abs(vrAt), abs(vrWh), vlZ0);
%     vrCorrA_VHA = calcCorr(abs(vrAt), abs(vrWh), vlZ0);
%     vrVafV_I = calcVAF(abs(vrV), vrI, vlZ0);
%     vrCorrV_I = calcCorr(abs(vrV), vrI, vlZ0);
%     vrVafVHA_I = calcVAF(abs(vrWh), vrI, vlZ0);
%     vrCorrVHA_I = calcCorr(abs(vrWh), vrI, vlZ0);
%     vrVafVA_I = calcVAF(abs(vrWt), vrI, vlZ0);
%     vrCorrVA_I = calcCorr(abs(vrWt), vrI, vlZ0);
%     vrVafA_I = calcVAF(abs(vrAt), vrI, vlZ0);
%     vrCorrAt_I = calcCorr(abs(vrAt), vrI, vlZ0);
%     
%     vrCorrAhf_I = calcCorr(abs(vrAhf), vrI1, vlZ0);
%     vrCorrAhv_I = calcCorr(abs(vrAhv), vrI1, vlZ0);
%     vrCorrAvf_I = calcCorr(abs(vrAvf), vrI1, vlZ0);
%     vrCorrAt_Ahv = calcCorr(abs(vrAt), abs(vrAhv1), vlZ0);
%     vrCorrAhm_I = calcCorr(abs(vrAhm), vrI1, vlZ0);
% 
%     vrVafI_D = calcVAF(vrI, vrD, vlZ0);
%     vrCorrI_D = calcCorr(vrI, vrD, vlZ0);    
%     vrVafV_D = calcVAF(abs(vrV), vrD, vlZ0);
%     vrCorrV_D = calcCorr(abs(vrV), vrD, vlZ0);    
%     [vrVafIV_D, Yf] = calcVAF([vrI1, vrV1, vrI1.*vrV1], vrD, vlZ0);
%     vrCorrIV_D = calcCorr(Yf, vrD, vlZ0);    
%     
%     vrVafI_DHA = calcVAF(vrI, abs(vrDAh), vlZ0);
%     vrCorrI_DHA = calcCorr(vrI, abs(vrDAh), vlZ0);    
%     vrVafVHA_DHA = calcVAF(abs(vrWh), abs(vrDAh), vlZ0);
%     vrCorrVHA_DHA = calcCorr(abs(vrWh), abs(vrDAh), vlZ0);    
%     [vrVafIVHA_DHA, Yf] = calcVAF([vrI1, vrWh1, vrI1.*vrWh1], vrDAh, vlZ0);
%     vrCorrIVHA_DHA = calcCorr(Yf, vrDAh1); 
%     
%     vrVafI_DA = calcVAF(vrI1, vrDA1);
%     vrCorrI_DA = calcCorr(vrI1, vrDA1);    
%     vrVafVA_DA = calcVAF(vrWt1, vrDA1);
%     vrCorrVA_DA = calcCorr(vrWt1, vrDA1);    
%     [vrVafIVA_DA, Yf] = calcVAF([vrI1, vrWt1, vrI1.*vrWt1], vrDA1);
%     vrCorrIVA_DA = calcCorr(Yf, vrDA1);     
%     
%     vrCorrV_A = calcCorr(vrV1, vrAt1);
%     vrCorrVHA_A = calcCorr(vrWh1, vrAt1);
%     vrCorrD_A = calcCorr(vrD1, vrAt1);
%     vrVafVHA_A = calcVAF(vrWh1, vrAt1);
%     vrVafV_A = calcVAF(vrV1, vrAt1);
%         
%     S = appendField(S, vrVafA_VHA, vrCorrA_VHA, vrVafV_I, vrCorrV_I, ...
%             vrVafVHA_I, vrCorrVHA_I, vrVafVA_I, vrCorrVA_I, vrVafA_I, ...
%             vrCorrAt_I, vrCorrAhf_I, vrCorrAhv_I, vrCorrAvf_I, ...
%             vrCorrAt_Ahv, vrCorrAhm_I, vrVafI_D, vrCorrI_D, vrVafV_D, ...
%             vrCorrV_D, vrVafIV_D, vrCorrIV_D, vrVafI_DHA, ...
%             vrCorrI_DHA, vrVafVHA_DHA, vrCorrVHA_DHA, vrVafIVHA_DHA, ...
%             vrCorrIVHA_DHA, vrVafI_DA, vrCorrI_DA, vrVafVA_DA, ...   
%             vrCorrVA_DA, vrVafIVA_DA, vrCorrIVA_DA, vrCorrV_A, ...
%             vrCorrVHA_A, vrCorrD_A, vrVafVHA_A, vrVafV_A, vrCorrV, vrCorrD, vrCorrI);            
end
S.vlZone = logical(S.vlZone);
end