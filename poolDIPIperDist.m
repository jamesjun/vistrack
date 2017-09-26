function S = poolDIPIperDist(vsTrial, iAnimal, mlMask)
%Distance to landmark and IPI

pixpercm = 1053.28/(sqrt(2)*100);
%landmark locations
xyf = [789, 681];    rf = 1; %cm, radius
xy1 = [966, 418];    r1 = 2.216*2.54/2; %*1.1222; %cm, radius
xy2 = [975, 790];    r2 = 3.545*2.54/2; %*1.1222; %cm, radius
xy3 = [604, 799];    r3 = 4*2.54/2; %cm, radius
xy4 = [600, 428];    r4 = 3*2.54/2; %cm, radius

if nargin >= 2 && ~isempty(iAnimal)
    viAnimal = poolVecFromStruct(vsTrial, 'iAnimal');
    vsTrial = vsTrial(viAnimal == iAnimal);
end

calcD1 = @(x,y)dist2square((x-xy1(1))/pixpercm, (y-xy1(2))/pixpercm, r1);
calcD2 = @(x,y)dist2square((x-xy2(1))/pixpercm, (y-xy2(2))/pixpercm, r2);
calcD3 = @(x,y)sqrt((x-xy3(1)).^2 + (y-xy3(2)).^2) / pixpercm - r3;        
calcD4 = @(x,y)sqrt((x-xy4(1)).^2 + (y-xy4(2)).^2) / pixpercm - r4;        
calcDf = @(x,y)sqrt((x-xyf(1)).^2 + (y-xyf(2)).^2) / pixpercm - rf;        

S.vrX = [];
S.vrY = [];
S.vrV = [];
S.vrI = [];
S.vrD = [];
S.vrDI = [];
S.vrD1 = []; %dist to LM1 (cm)
S.vrD2 = []; %dist to LM2 (cm)
S.vrD3 = []; %dist to LM3 (cm)
S.vrD4 = []; %dist to LM4 (cm)
S.vrDf = []; %dist to Food (cm)
S.vrCorrV = [];
S.vrCorrD = [];
S.vrCorrI = [];
S.img0 = vsTrial(1).img0;
S.xy0 = vsTrial(1).xy0; 

for iTrial=1:numel(vsTrial)
    Strial = vsTrial(iTrial);
    
    vtEOD = poolVecFromStruct(Strial, 'vtEOD');
    vrXr = poolVecFromStruct(Strial, 'vrX');
    vrYr = poolVecFromStruct(Strial, 'vrY');
    vrTr = poolVecFromStruct(Strial, 'TEOD');
    vrVr = poolVecFromStruct(Strial, 'VEL');

    % interpolate the location of the EOD pulse
    vrX = interp1(vrTr, vrXr, vtEOD, 'spline');
    vrY = interp1(vrTr, vrYr, vtEOD, 'spline');    
    vrV = interp1(vrTr, vrVr, vtEOD, 'spline');
    vrD = sqrt(differentiate3(vrX).^2 + differentiate3(vrY).^2) / pixpercm;
    vrI = differentiate3(vtEOD);
    vrDI = differentiate3(vrI);
    vrD1 = calcD1(vrX, vrY);
    vrD2 = calcD2(vrX, vrY);
    vrD3 = calcD3(vrX, vrY);
    vrD4 = calcD4(vrX, vrY);
    vrDf = calcDf(vrX, vrY);
    
    if nargin >= 3
        vlRegion = mlMask(sub2ind(size(mlMask), round(vrY), round(vrX)));
        vrX = vrX(vlRegion);
        vrY = vrY(vlRegion);
        vrI = vrI(vlRegion);
        vrD = vrD(vlRegion);
        vrDI = vrDI(vlRegion);
        vrD1 = vrD1(vlRegion);
        vrD2 = vrD2(vlRegion);
        vrD3 = vrD3(vlRegion);
        vrD4 = vrD4(vlRegion);
        vrDf = vrDf(vlRegion);
        vrV = vrV(vlRegion);
    end

    S.vrX = [S.vrX; vrX(:)];
    S.vrY = [S.vrY; vrY(:)];
    S.vrV = [S.vrV; vrV(:)];
    S.vrI = [S.vrI; vrI(:)];
    S.vrD = [S.vrD; vrD(:)];
    S.vrDI = [S.vrDI; vrDI(:)];
    S.vrD1 = [S.vrD1; vrD1(:)];
    S.vrD2 = [S.vrD2; vrD2(:)];
    S.vrD3 = [S.vrD3; vrD3(:)];
    S.vrD4 = [S.vrD4; vrD4(:)];
    S.vrDf = [S.vrDf; vrDf(:)];
    S.vrCorrV(end+1) = calcCorrTau(vrV);
    S.vrCorrD(end+1) = calcCorrTau(vrD);
    S.vrCorrI(end+1) = calcCorrTau(vrI);
end

