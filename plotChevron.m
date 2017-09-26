function h = plotChevron(XI, YI, vrColor, ANGLE, scale)
% XI: [x_tip, x_tail], YI: [y_tip, y_tail]

if nargin < 3
    vrColor = [1 0 0]; %plot red by default
end
if nargin < 4
    ANGLE = 60; %chevron angle in degrees
end
if nargin < 5
    scale = 1;
end

% tip of the chevron
vrX(2) = XI(1);
vrY(2) = YI(1);

% Rotate vectors
vec0 = [XI(2) - XI(1), YI(2) - YI(1)];
vecP = rotatexy(vec0, ANGLE/2) * scale;
vecN = rotatexy(vec0, -ANGLE/2) * scale;
vrX(1) = vecP(1) + vrX(2);
vrX(3) = vecN(1) + vrX(2);
vrY(1) = vecP(2) + vrY(2);
vrY(3) = vecN(2) + vrY(2);

% plot chevron
h = plot(vrX, vrY, '-', 'color', vrColor);

end

function [ xyp ] = rotatexy( xy, ang )
%ROTATEXY rotate a vector wrt origin, ang in degree

xy = xy(:);
CosA = cos(deg2rad(ang));
SinA = sin(deg2rad(ang));

M = [CosA, -SinA; SinA, CosA];
xyp = M * xy;

end