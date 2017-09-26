function [Xr, Yr] = rotXY(X, Y, xy0, angRot, dimm)
% rotate x,y coordinates by the angle, given in degrees
rotMat = rotz(angRot);    rotMat = rotMat(1:2, 1:2);

if nargin < 5
    dimm = 1;
end

if min(size(X)) > 1
    if dimm == 2
        X = X';
        Y = Y';
    end
    Xr = zeros(size(X));
    Yr = zeros(size(Y));
    for iCol=1:size(X,2)
        [Xr(:,iCol), Yr(:,iCol)] = ...
            rotXY(X(:,iCol), Y(:,iCol), xy0, angRot);
    end
    if dimm == 2
        Xr = Xr';
        Yr = Yr';
    end    
    return;
end
    
    
mrXY = [X(:) - xy0(1), Y(:) - xy0(2)] * rotMat;
Xr = mrXY(:,1) + xy0(1);
Yr = mrXY(:,2) + xy0(2);
