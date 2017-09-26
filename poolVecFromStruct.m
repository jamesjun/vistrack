function vr = poolVecFromStruct(vs, fieldname, mlMask, iAnimal)

% filter by animals
if nargin >= 4
    if ~isempty(iAnimal) && iAnimal > 0
        viAnimal = poolVecFromStruct(vs, 'iAnimal');
        vs = vs(viAnimal == iAnimal);
    end
end

vr = [];
for i=1:numel(vs)
    if iscell(vs)
        vs1 = vs{i};
    else
        vs1 = vs(i);
    end
    vr1 = getfield(vs1, fieldname);

    if nargin >= 3 && ~isempty(mlMask)
        vr1 = vr1(getRegion(vs1, mlMask));
    end

    vr = [vr; vr1(:)];
end
end