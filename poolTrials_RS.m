function S = poolTrials_RS(vsTrial, iAnimal, lim, iZone)
%Distance to landmark and IPI

if nargin < 3, lim = []; end %no trim
if nargin < 4, iZone = 1; end %active zone

pixpercm = 1053.28/(sqrt(2)*100);
%landmark locations
xy0 = vsTrial(1).xy0;
xyf = [789, 681];    rf = 1; %cm, radius
xy1 = [966, 418];    r1 = 2.216*2.54/2; %*1.1222; %cm, radius
xy2 = [975, 790];    r2 = 3.545*2.54/2; %*1.1222; %cm, radius
xy3 = [604, 799];    r3 = 4*2.54/2; %cm, radius
xy4 = [600, 428];    r4 = 3*2.54/2; %cm, radius
calcSpeed = @(x,y)sqrt(differentiate5(x).^2 + differentiate5(y).^2) * 100 / pixpercm;
calcAng = @(x,y)cart2pol(differentiate5(x), differentiate5(y)); %tangent angle in rad

%Select animal-specific trials
if nargin >= 2
    if ~isempty(iAnimal)
        if iAnimal > 0
            viAnimal = poolVecFromStruct(vsTrial, 'iAnimal');
            vsTrial = vsTrial(viAnimal == iAnimal);
        end
    end
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
        'vrX', 'vrY', 'VEL', 'EODR', 'ANG', 'HANG', 'AVEL', 'HAVEL', 'mrX', 'mrY');
    
    vrX = Strial.vrX;
    vrY = Strial.vrY;
    vrXm = Strial.mrX(:,4); %midpoint
    vrYm = Strial.mrY(:,4);
    vrSh = calcSpeed(vrX, vrY);
    vrSm = calcSpeed(vrXm, vrYm);
    vrSr = vrSh ./ vrSm; %speed ratio
    vrAv = calcAng(vrX, vrY)';
    vrAm = calcAng(vrXm, vrYm)';
    vrAf = cart2pol(xyf(1)-vrX, xyf(2)-vrY);
    vrV = Strial.VEL' / pixpercm;
    vrR = Strial.EODR;
    vrAt = Strial.ANG;
    vrAh = Strial.HANG;
    vrWt = Strial.AVEL;
    vrWh = Strial.HAVEL;
    vrD0 = calcD0(vrX, vrY);
    vrD1 = calcD1(vrX, vrY);
    vrD2 = calcD2(vrX, vrY);
    vrD3 = calcD3(vrX, vrY);
    vrD4 = calcD4(vrX, vrY);
    vrDf = calcDf(vrX, vrY);    
    vlZone = isZone(vrX, vrY, Strial.xy0);
    vpBack = mean(vrV(vlZone)<0);
    vrAe = calcAngErr(vrAf, vrAh); 
    
    S1 = makeStruct(vrD0, vrD1, vrD2, vrD3, vrD4, vrDf, vlZone);
    vlZ0 = getZone(S1, iZone);
    vrCorrV_At = calcCorr(abs(vrV), abs(vrAt), vlZ0);
    vrCorrVh_At = calcCorr(abs(vrSh), abs(vrAt), vlZ0);
    vrCorrAe_At = calcCorr(abs(vrAe), abs(vrAt), vlZ0);

    S = appendField(S, vrX, vrY, vrV, vrR, vrSr, vrSm, vrSh, ...
            vrAt, vrAh, vrAv, vrAf, vrWt, vrWh, vrAm, vrAe, ...
            vrD0, vrD1, vrD2, vrD3, vrD4, vrDf, vlZone, vpBack, ...
            vrCorrV_At, vrCorrVh_At, vrCorrAe_At, vlZ0);
end
S.vlZone = logical(S.vlZone);
end