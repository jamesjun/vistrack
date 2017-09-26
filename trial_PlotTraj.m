function trial_PlotTraj(S, varargin)
xy_names = {'CoM', 'head', 'head-mid', 'mid', 'tail-mid', 'tail'};

P = struct('vrColor', [0 0 1], 'timeInterval', 50, 'radLim', [1 10], ...
    'fFadePath', 1, 'fFadeTime', 1, 'fAccum', 0, 'nxy', [] ,...
    'angRot', -1.1590, 'interpStep', .2, 'img0', []);
if nargin > 1
    P = appendStruct(P, varargin{:});
end

if ~isempty(P.angRot)
    [mrX, mrY] = trial_correctRot(S, P.angRot);
else
    mrX = S.mrX;
    mrY = S.mrY;
    P.angRot = 0;
end
if isempty(P.img0)
    img0 = S.img0;
else
    img0 = P.img0;
end

if ~P.fAccum, 
    img0 = imrotate(img0, P.angRot);
    imshow(imadjust(img0)); 
end
hold on;
nTime = size(mrX, 1);
viTime = 1:nTime;
vrFade = interp1(viTime([1 end]), [0 1], viTime);
viTime = nTime:-P.timeInterval:1;
viTime = fliplr(viTime);

if isempty(P.nxy), nxy = size(mrX,2);
else nxy = P.nxy; end
nPath = numel(2:P.interpStep:nxy);
vrRadius = interp1([1 nPath], P.radLim, 1:nPath);
vrFadePath = interp1([1 nPath], [0 1], 1:nPath);


for iTime=viTime
    vrXp = fliplr(interp1(2:nxy, mrX(iTime, 2:nxy), 2:P.interpStep:nxy, 'spline'));
    vrYp = fliplr(interp1(2:nxy, mrY(iTime, 2:nxy), 2:P.interpStep:nxy, 'spline'));
    if P.fFadeTime
        vrColor = vrFade(iTime) * P.vrColor + (1-vrFade(iTime)) * [1 1 1];
    else
        vrColor = P.vrColor;
    end
    for iPath = 1:nPath
        if P.fFadePath
            vrColor1 = vrFadePath(iPath) * vrColor + (1-vrFadePath(iPath)) * [1 1 1];
        else
            vrColor1 = vrColor;
        end
        plot(vrXp(iPath), vrYp(iPath), 'o', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', vrColor1, 'MarkerSize', vrRadius(iPath)); %head
    end
end
% plot(mrX(:,2), mrY(:,2), 'r.'); hold on;
% plot(mrX(:,2), mrY(:,2), 'b-'); hold on;
%plot(mrX(:,2), mrY(:,2), 'k:'); hold on; %midpoint
title(S.dataID);