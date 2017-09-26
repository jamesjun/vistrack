function vrD = dist2square(vrX, vrY, r)
% square centre located at (0,0) and edges parallel to the x,y axis
% r is the half-length of the square
% calculate the distance to the nearest edge

% [th, rho] = cart2pol(x,y); %rad
calcD = @(x,y,vl)sqrt((vrX(vl) - x).^2 + (vrY(vl) - y).^2);     

vlX1 = vrX < -r;
vlX3 = vrX >= r;
vlX2 = ~vlX1 & ~vlX3;

vlY1 = vrY < -r;
vlY3 = vrY >= r;
vlY2 = ~vlY1 & ~vlY3;

vrD = zeros(size(vrX));

vrD(vlX1 & vlY1) = calcD(-r,-r, vlX1 & vlY1);
vrD(vlX1 & vlY3) = calcD(-r,+r, vlX1 & vlY3);
vrD(vlX3 & vlY1) = calcD(+r,-r, vlX3 & vlY1);
vrD(vlX3 & vlY3) = calcD(+r,+r, vlX3 & vlY3);

vrD(vlX2 & vlY3) = vrY(vlX2 & vlY3) - r;
vrD(vlX2 & vlY1) = -r - vrY(vlX2 & vlY1);
vrD(vlX1 & vlY2) = -r - vrX(vlX1 & vlY2);
vrD(vlX3 & vlY2) = vrX(vlX3 & vlY2) - r;


