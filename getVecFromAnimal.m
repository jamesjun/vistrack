function vrVec = getVecFromAnimal(vsTrial, strVec, iAnimal, mlMask)

if ~isempty(iAnimal)    
    viAnimal = poolVecFromStruct(vsTrial, 'iAnimal');
    vsTrial = vsTrial(viAnimal == iAnimal);
end

vrVec = poolVecFromStruct(vsTrial, strVec);

if nargin >= 4
    vrVec = vrVec(getRegion(vsTrial, mlMask));
end