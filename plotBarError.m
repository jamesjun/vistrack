function plotBarError(mrMean, mrSem, csColor, XTickLabel)
% mrMean: n by m, group sof m, and n categories
nX = size(mrMean,1);
h = bar(mrMean, .75, 'EdgeColor', 'none'); hold on;

% assign bar colors
if nargin < 2
    mrSem = [];
end
if nargin < 3
    csColor = [];
end
if nargin < 4
    XTickLabel = [];
end

if ~isempty(mrSem);    
    for ih=1:numel(h)
        vrX = mean(get(get(h(ih), 'Children'), 'XData'), 1);
        colorErr = 'r.'; %'.', 'color', csColor{ih}
        errorbar(vrX, mrMean(:,ih), mrSem(:,ih), colorErr);        
    end
end


if ~isempty(csColor)
    for ih=1:numel(h)
        set(h(ih), 'FaceColor', csColor{ih});
    end
end


if ~isempty(XTickLabel)
    % assign xlabels
    set(gca, 'XLim', [.5, nX+.5]);
    set(gca, 'XTick', 1:nX);
    if nargin >= 4
        set(gca, 'XTickLabel', XTickLabel);
    end    
end

end %func