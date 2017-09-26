function [mrDx, mrDy, csDist] = calcDist2Fish(S)
% pass a trial object and get the distance matrix to fish

% do not fix this location. for random moving code the object location in
% the GUI tool
xy0 = vsTrial(1).xy0;
xyf = [789, 681];    rf = 1; %cm, radius
xy1 = [966, 418];    r1 = 2.216*2.54/2; %*1.1222; %cm, radius
xy2 = [975, 790];    r2 = 3.545*2.54/2; %*1.1222; %cm, radius
xy3 = [604, 799];    r3 = 4*2.54/2; %cm, radius
xy4 = [600, 428];    r4 = 3*2.54/2; %cm, radius


calcD0 = @(x,y)sqrt((x-xy0(1)).^2 + (y-xy0(2)).^2) / pixpercm;        
calcD1 = @(x,y)dist2square((x-xy1(1))/pixpercm, (y-xy1(2))/pixpercm, r1);
calcD2 = @(x,y)dist2square((x-xy2(1))/pixpercm, (y-xy2(2))/pixpercm, r2);
calcD3 = @(x,y)sqrt((x-xy3(1)).^2 + (y-xy3(2)).^2) / pixpercm - r3;        
calcD4 = @(x,y)sqrt((x-xy4(1)).^2 + (y-xy4(2)).^2) / pixpercm - r4;        
calcDf = @(x,y)sqrt((x-xyf(1)).^2 + (y-xyf(2)).^2) / pixpercm - rf;

csDist = {'Center', 'LM1', 'LM2', 'LM3', 'LM4'};