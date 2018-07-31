function [VISITCNT, TIMECNT] = calcVisitDensity(I0, TC, XC, YC, TRAJ_NFILT, evalstr)
R = 5; %2 4 8 10 16 20* 32
timeout = 25; %25ms timeout
step = 1;

if exist('evalstr', 'var')
    eval(evalstr);
end

% interpolate X,Y
TCi = TC(1):(1/100):TC(end);
Xi = interp1(TC, filtPos(XC, TRAJ_NFILT), TCi, 'spline', 'extrap'); 
Yi = interp1(TC, filtPos(YC, TRAJ_NFILT), TCi, 'spline', 'extrap'); 

% tic;
XLIM = [1, size(I0,2)];     YLIM = [1, size(I0,1)];

% build cache
[XX, YY] = meshgrid(-R:1:R);
XX = XX(:);                 YY = YY(:);
RR = XX .^ 2 + YY .^ 2;
XX(RR > R^2) = [];          YY(RR > R^2) = [];

%count time spent
TIMECNT = zeros(size(I0)); %1 means 1msec stay

%count time spent
VISITCNT = zeros(size(I0)); %1 means 1msec stay
HISTORY = zeros(size(I0), 'uint8');
for i=1:step:numel(Xi)        
    % Get nearby indices
    X = XX + round(Xi(i));      Y = YY + round(Yi(i));
    IDXKILL = X < XLIM(1) | X > XLIM(2) | Y < YLIM(1) | Y > YLIM(2);
    X(IDXKILL) = [];            Y(IDXKILL) = [];
    
    TIMECNT(Y,X) = TIMECNT(Y,X) + step;        
    idx = sub2ind(size(I0), Y, X);    
    idx2 = idx(HISTORY(idx) == 0); % keep non visited index
    if ~isempty(idx2)
        VISITCNT(idx2) = VISITCNT(idx2) + 1;
    end    
    HISTORY(idx) = timeout;
    HISTORY=HISTORY-step;
end


% make grid 2 homogenize
VISITCNT(2:2:end, 2:2:end) = VISITCNT(1:2:end, 1:2:end);
VISITCNT(1:2:end, 2:2:end) = VISITCNT(1:2:end, 1:2:end);
VISITCNT(2:2:end, 1:2:end) = VISITCNT(1:2:end, 1:2:end);

if nargout == 0
    figure; 
    subplot 121; 
    imshow(rgbmix(I0, imgray2rgb((TIMECNT))));   title('Time density map');
    
    subplot 122; 
    imshow(rgbmix(I0, imgray2rgb((VISITCNT))));  title('Visit density map');    
end
end %func