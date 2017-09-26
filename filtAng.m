function [yw, y] = filtAng( x, n )
%x is radian
%FILTPOS Apply smoothing filter to the orientation angle
%James J. Jun. Malerlab. 2012 Nov
%
%[Variables]
%yw: wrapped orientation angles (ranges between -180 to 180 deg)
%y: unwrapped orientation angles

if nargin < 2, n = 15; end

% x = rad2deg(unwrap(deg2rad(x), pi));  
x = unwrap(x, pi);  

n1 = round(n/2); if mod(n1,2), n1 = n1 + 1; end
y1 = medfilt1(x,n1);
y = filtfilt(ones(1,n), n, y1);
y(1:n) = median(x(1:n)); %fix the start

yw = mod(y, 2*pi); 
yw(yw>pi) = yw(yw>pi) - 2*pi;
end

