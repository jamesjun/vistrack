function S = poolDistPerIPI(vsTrial, iAnimal, mlMask)
pixpercm = 1053.28/(sqrt(2)*100);

if nargin >= 2
    if ~isempty(iAnimal) && iAnimal > 0
        viAnimal = poolVecFromStruct(vsTrial, 'iAnimal');
        vsTrial = vsTrial(viAnimal == iAnimal);
    end
end

midVal = @(x)(x(1:end-1) + x(2:end))/2;

S.vrD = [];
S.vrX = [];
S.vrY = [];
S.vrI = [];
S.vrV = [];
S.vrDI = [];
S.vrDIn = [];
S.vrDIs = [];
S.vrDIns = [];

for iTrial=1:numel(vsTrial)
    Strial = vsTrial(iTrial);
    
    vtEOD = poolVecFromStruct(Strial, 'vtEOD');
    vrXr = poolVecFromStruct(Strial, 'vrX');
    vrYr = poolVecFromStruct(Strial, 'vrY');
    vrTr = poolVecFromStruct(Strial, 'TEOD');
    vrVr = poolVecFromStruct(Strial, 'VEL');

    % interpolate the location of the EOD pulse
    vrXe = interp1(vrTr, vrXr, vtEOD, 'spline');
    vrYe = interp1(vrTr, vrYr, vtEOD, 'spline');
    vrD = sqrt(diff(vrXe).^2 + diff(vrYe).^2) / pixpercm;        
    
    vrX = midVal(vrXe);
    vrY = midVal(vrYe);
    vrI = diff(vtEOD);
    vrDI = diff(vrI);
    vrDI(end+1) = vrDI(end);
    vrDIn = vrDI ./ vrI;
    
    vrT = midVal(vtEOD);
    vrV = interp1(vrTr, vrVr, vrT, 'spline');    
    
    if nargin >= 3
        vlRegion = mlMask(sub2ind(size(mlMask), round(vrY), round(vrX)));
        vrD = vrD(vlRegion);
        vrX = vrX(vlRegion);
        vrY = vrY(vlRegion);
        vrI = vrI(vlRegion);
        vrV = vrV(vlRegion);
        vrDI = vrDI(vlRegion);
        vrDIn = vrDIn(vlRegion);
    end

    
    S.vrD = [S.vrD; vrD(:)];
    S.vrX = [S.vrX; vrX(:)];
    S.vrY = [S.vrY; vrY(:)];
    S.vrI = [S.vrI; vrI(:)];
    S.vrDI = [S.vrDI; vrDI(:)];
    S.vrV = [S.vrV; vrV(:)];
    S.vrDIn = [S.vrDIn; vrDIn];
    S.vrDIs = [S.vrDIs(:); calcDistAsym(vrDI)];
    S.vrDIns = [S.vrDIns(:); calcDistAsym(vrDIn)];
end