function mrHist = trial_PlotTrajBitmap(S, varargin)
xy_names = {'CoM', 'head', 'head-mid', 'mid', 'tail-mid', 'tail'};

P = struct('vrColor', [0 0 1], 'timeInterval', 50, 'radLim', [10 1], ...
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
height = size(img0, 1);
width = size(img0, 2);
mrHist = zeros(size(img0));
for iTime=viTime
    vrXp = interp1(2:nxy, mrX(iTime, 2:nxy), 2:P.interpStep:nxy, 'spline');
    vrYp = interp1(2:nxy, mrY(iTime, 2:nxy), 2:P.interpStep:nxy, 'spline');
    mlMask = false(size(img0));    
    for iPath = 1:nPath
        [vrXc, vrYc] = getCircle([vrXp(iPath), vrYp(iPath)], vrRadius(iPath));
        mlMask = mlMask | poly2mask(vrXc, vrYc, height, width);
    end
    % accumulate mask
    mrHist(mlMask) + mrHist(mlMask) + 1;
end
if nargout == 0
    figure; imagesc(mrHist);
    title(S.dataID);
end