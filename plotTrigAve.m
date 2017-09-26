function [Ymean, Ystd, YMAT, Ypre] = plotTrigAve(Y, IDXTRIG, IDXLIM, IDXPRELIM)
%[Ymean, Ystd, YPLOT, xinterp] = plotTriggeredAveraging(T, Y, Ttrig, twin, tdur, baselineRange)
%tdur: stimulus duration
% averaging onset time is proven to be accurate

if nargin < 3    
    IDXLIM = [-500, 500];
end
if nargin < 4
    IDXPRELIM = [];
end

%initialize
NWIN = IDXLIM(2) - IDXLIM(1) + 1;
nTriggers = numel(IDXTRIG);
YMAT = zeros(nTriggers, NWIN);
Ypre = zeros(nTriggers, 1);

%construct window
viKill = [];
for iTrig=1:nTriggers
    try
        idx0 = IDXTRIG(iTrig);
        IDX = idx0 + [IDXLIM(1):IDXLIM(2)];
        if ~isempty(IDXPRELIM)
            IDXPRE = idx0 + [IDXPRELIM(1):IDXPRELIM(2)];
            Ypre(iTrig) = mean(Y(IDXPRE));
            YMAT(iTrig, :) = Y(IDX) - Ypre(iTrig);
        else
            YMAT(iTrig, :) = Y(IDX);
        end
    catch
        viKill(end+1) = iTrig;
    end
end
YMAT(viKill,:) = [];

% YMAT = abs(YMAT);

Ymean = mean(YMAT);
Ystd = std(YMAT);

% show all
if nargout == 0
    figure;
    
    XPLOT = (IDXLIM(1):IDXLIM(2));    
    plot(XPLOT, YMAT', '-', 'Color', .85*ones(3,1)); hold on;
    plot(XPLOT, Ymean, 'b-');    
    plot(XPLOT, Ymean+Ystd, 'b:');    
    plot(XPLOT, Ymean-Ystd, 'b:');    
    xlabel('IDX triggered at 0');

    axis([IDXLIM(1), IDXLIM(2), min(YMAT(:)), max(YMAT(:))]);
    plot([0 0], get(gca, 'YLim'), 'r-');
     grid on;
end
