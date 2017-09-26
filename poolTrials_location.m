function S = poolTrials_location(vsTrial, iAnimal)
%Distance to landmark and IPI

pixpercm = 1053.28/(sqrt(2)*100);
%landmark locations
xy0 = vsTrial(1).xy0;
xyf = [789, 681];    rf = 1; %cm, radius
xy1 = [966, 418];    r1 = 2.216*2.54/2; %*1.1222; %cm, radius
xy2 = [975, 790];    r2 = 3.545*2.54/2; %*1.1222; %cm, radius
xy3 = [604, 799];    r3 = 4*2.54/2; %cm, radius
xy4 = [600, 428];    r4 = 3*2.54/2; %cm, radius

if nargin >= 2 && ~isempty(iAnimal)
    viAnimal = poolVecFromStruct(vsTrial, 'iAnimal');
    vsTrial = vsTrial(viAnimal == iAnimal);
end

calcD0 = @(x,y)sqrt((x-xy0(1)).^2 + (y-xy0(2)).^2) / pixpercm;        
calcD1 = @(x,y)dist2square((x-xy1(1))/pixpercm, (y-xy1(2))/pixpercm, r1);
calcD2 = @(x,y)dist2square((x-xy2(1))/pixpercm, (y-xy2(2))/pixpercm, r2);
calcD3 = @(x,y)sqrt((x-xy3(1)).^2 + (y-xy3(2)).^2) / pixpercm - r3;        
calcD4 = @(x,y)sqrt((x-xy4(1)).^2 + (y-xy4(2)).^2) / pixpercm - r4;        
calcDf = @(x,y)sqrt((x-xyf(1)).^2 + (y-xyf(2)).^2) / pixpercm - rf;        


S.vrX = [];
S.vrY = [];
S.vrV = [];
S.vrR = [];
S.vrA = [];
S.vrAV = [];
S.vrD0 = []; %dist from centre
S.vrD1 = []; %dist to LM1 (cm)
S.vrD2 = []; %dist to LM2 (cm)
S.vrD3 = []; %dist to LM3 (cm)
S.vrD4 = []; %dist to LM4 (cm)
S.vrDf = []; %dist to Food (cm)
S.vlZone = []; %active zone
S.img0 = vsTrial(1).img0;
S.xy0 = vsTrial(1).xy0; 

for iTrial=1:numel(vsTrial)
    Strial = vsTrial(iTrial);
    
    vrX = poolVecFromStruct(Strial, 'vrX');
    vrY = poolVecFromStruct(Strial, 'vrY');
    vrV = poolVecFromStruct(Strial, 'VEL');
    vrR = poolVecFromStruct(Strial, 'EODR');
    vrA = poolVecFromStruct(Strial, 'ANG');
    vrAV = poolVecFromStruct(Strial, 'AVEL');
    vrD0 = calcD0(vrX, vrY);
    vrD1 = calcD1(vrX, vrY);
    vrD2 = calcD2(vrX, vrY);
    vrD3 = calcD3(vrX, vrY);
    vrD4 = calcD4(vrX, vrY);
    vrDf = calcDf(vrX, vrY);    
    vlZone = isZone(Strial);

    S.vrX = [S.vrX; vrX(:)];
    S.vrY = [S.vrY; vrY(:)];
    S.vrV = [S.vrV; vrV(:)];
    S.vrR = [S.vrR; vrR(:)];
    S.vrA = [S.vrA; vrA(:)];
    S.vrAV = [S.vrAV; vrAV(:)];
    S.vrD0 = [S.vrD0; vrD0(:)];
    S.vrD1 = [S.vrD1; vrD1(:)];
    S.vrD2 = [S.vrD2; vrD2(:)];
    S.vrD3 = [S.vrD3; vrD3(:)];
    S.vrD4 = [S.vrD4; vrD4(:)];
    S.vrDf = [S.vrDf; vrDf(:)];
    S.vlZone = [S.vlZone; vlZone(:)];
end
S.vlZone = logical(S.vlZone);
end

function vl = isZone(S)
angRot = -1.1590; %deg
rectCrop = [493 1083 312 902];

% rotational correction
rotMat = rotz(angRot);    rotMat = rotMat(1:2, 1:2);
mrXY = [S.vrX(:) - S.xy0(1), S.vrY(:) - S.xy0(2)] * rotMat;
vrX = mrXY(:,1) + S.xy0(1);
vrY = mrXY(:,2) + S.xy0(2);

vl = vrX >= rectCrop(1) & vrX < rectCrop(2) ...
   & vrY >= rectCrop(3) & vrY < rectCrop(4);

end

