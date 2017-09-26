function [xyStart, xyFood, angRot] = getCoordinatesForAnimal(fishID)

switch fishID
    case 'A'
        xyStart = [55, 50]; xyFood = [0 -10]; angRot = 0;
    case 'B'
        xyStart = [50, -55]; xyFood = [-10 0]; angRot = 90;
    case 'C'
        xyStart = [-55, -50]; xyFood = [0 10]; angRot = 180;
    case 'D'
        xyStart = [-50, 55]; xyFood = [10 0]; angRot = 270;
end