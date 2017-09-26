function [vl, rectCrop] = isZone(vrX, vrY, xy0)
angRot = -1.1590; %deg
rectCrop0 = [493 1083 312 902];
dxy = xy0 - [787.0169, 605.6858];
rectCrop = rectCrop0 + [dxy(1), dxy(1), dxy(2), dxy(2)];

% rotational correction
rotMat = rotz(angRot);    rotMat = rotMat(1:2, 1:2);
mrXY = [vrX(:) - xy0(1), vrY(:) - xy0(2)] * rotMat;
vrX = mrXY(:,1) + xy0(1);
vrY = mrXY(:,2) + xy0(2);

vl = vrX >= rectCrop(1) & vrX < rectCrop(2) ...
   & vrY >= rectCrop(3) & vrY < rectCrop(4);

end