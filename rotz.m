function [ mrRot ] = rotz( angdeg )
%JJJ implemented
%angdeg: angle in degrees
angrad = angdeg / 180 * pi;

mrRot = [cos(angrad), -sin(angrad), 0; ...
         sin(angrad), cos(angrad), 0; ...
        0, 0, 1];

end

