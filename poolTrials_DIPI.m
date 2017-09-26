function S = poolTrials_DIPI(vsTrial, iAnimal, lim)
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

if nargin < 2
    iAnimal = []; %time index
end
if nargin < 3
    lim = []; %time index
end

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
    vrX = interp1(vrTr, Strial.vrX, vtEOD, 'spline');
    vrY = interp1(vrTr, Strial.vrY, vtEOD, 'spline');    
    vrV = interp1(vrTr, Strial.VEL, vtEOD, 'spline') / pixpercm;
    vrAcc = interp1(vrTr, Strial.ACC, vtEOD, 'spline') / pixpercm;
    
    % interpolate the location of the EOD pulse
    vrI = differentiate3(vtEOD);
    vrDI = differentiate3(vrI);
    viT = findDIsac(vrDI);
    
    
    %angle
    vrAvr = unwrap(cart2pol(differentiate3(Strial.vrX), differentiate3(Strial.vrY))); %tangential angle
    vrWvr = differentiate5(vrAvr);
    vrAv = interp1(vrTr, vrAvr, vtEOD, 'spline');    
    vrAh = interp1(vrTr, unwrap(Strial.HANG), vtEOD, 'spline');
    vrAt = interp1(vrTr, Strial.ANG, vtEOD, 'spline');    
    vrWv = interp1(vrTr, vrWvr, vtEOD, 'spline');
    vrWh = interp1(vrTr, Strial.HAVEL, vtEOD, 'spline');
    vrWt = interp1(vrTr, Strial.AVEL, vtEOD, 'spline');
    vrALt = interp1(vrTr, differentiate5(Strial.AVEL), vtEOD, 'spline'); %ang accel
    vrAf = cart2pol(xyf(1)-vrX, xyf(2)-vrY);
    vrAhf = calcAngErr(vrAh, vrAf);
    vrAhv = calcAngErr(vrAh, vrAv);
    vrAvf = calcAngErr(vrAv, vrAf);
    
    %per IPI
    vrD = sqrt(differentiate3(vrX).^2 + differentiate3(vrY).^2) / pixpercm;  
    
    vrDAv = differentiate3(vrAv);
    vrDAh = differentiate3(vrAh);
    vrDAt = differentiate3(vrAt);    
    vrAv = wrap(vrAv);
    vrAh = wrap(vrAh);
%     vrDA = interp1(vrTr, poolVecFromStruct(Strial, 'AVEL'), vtEOD, 'spline');
    
    vrIz = zscore(vrI);
    vrIq = qscore(vrI);
    
    vrD0 = calcD0(vrX, vrY);
    vrD1 = calcD1(vrX, vrY);
    vrD2 = calcD2(vrX, vrY);
    vrD3 = calcD3(vrX, vrY);
    vrD4 = calcD4(vrX, vrY);
    vrDf = calcDf(vrX, vrY);
    vlZone = isZone(vrX, vrY, Strial.xy0);
    vrPmm = 1/mean(vrD(vlZone)*10);
    vrPmm_F15 = 1/mean(vrD(vrDf<15)*10);
    
    S = appendField(S, vrX, vrY, vrV, vrAcc, ...
            vrI, vrDI, vrIz, vrIq, ...
            vrD, vrDAt, vrDAh, vrDAv, vrPmm, vrPmm_F15, ...
            vrD0, vrD1, vrD2, vrD3, vrD4, vrDf, vlZone, ...
            vrAh, vrAt, vrAv, vrWh, vrWt, vrWv, vrALt, ...
            vrAf);
  
end
S.vlZone = logical(S.vlZone);
end