function mrY = plotBarErr(varargin)
%n by 3 (val, low, high)
nBoot = 100;

bootMean = @(x)[mean(x); bootci(nBoot, {@(y)mean(y), x}, 'type', 'per')];

n = nargin;
vrX = 1:n;

mrY = zeros(3, n);
for i=1:n
    if numel(varargin{i}) == 1
        mrY(:,i) = repmat(varargin{i}, [1,3]);
    else
        mrY(:,i) = bootMean(varargin{i});
    end
end
errorbar(vrX, mrY(1,:), mrY(1,:)-mrY(2,:), mrY(3,:)-mrY(1,:), 'r.'); hold on;

bar(vrX, mrY(1,:), .5);
set(gca, 'XLim', [.5 n+.5]);
set(gca, 'XTick', vrX);