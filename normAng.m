function vrA = normAng(vrA)
vrA = mod(vrA, 360);
vl = vrA > 180;
vrA(vl) = vrA(vl) - 360;

vl = vrA < -180;
vrA(vl) = vrA(vl) + 360;