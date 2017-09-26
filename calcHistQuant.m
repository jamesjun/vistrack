function [mlX_E, mlX_L] = calcHistQuant(mrX, viSession_E, viSession_L, quantLim, xi)
% select subset of data 

if nargin < 5
    xi = linspace(min(mrX(:)), max(mrX(:)), 32);
end

mrX_E = mrX(viSession_E,:);
mrX_L = mrX(viSession_L,:);
quant_E = quantile(mrX_E(:), quantLim);
quant_L = quantile(mrX_L(:), quantLim);
mlX_E = (mrX_E >= quant_E(1)) & (mrX_E < quant_E(2));
mlX_L = (mrX_L >= quant_L(1)) & (mrX_L < quant_L(2));

vrX_E = mrX_E(:);
vrX_L = mrX_L(:);

vnE = hist(vrX_E, xi); 
vnL = hist(vrX_L, xi);
vnE1 = hist(vrX_E(mlX_E), xi);
vnL1 = hist(vrX_L(mlX_L), xi);

strRange = sprintf('%0.1f~%0.1f%%', quantLim(1)*100, quantLim(2)*100);


% if nargout == 0
    hold on;
    h = bar(xi, [vnE; vnL]', 1, 'FaceColor', 'none');
    set(h(1), 'EdgeColor', 'r');    
    set(h(2), 'EdgeColor', 'b');
    ylabel('# Obs.');
%     legend({'Early all', 'Late all'}); xlabel('X'); ylabel('# Obs.');
    
    h = bar(xi, [vnE1; vnL1]', 1, 'EdgeColor', 'none');
    set(h(1), 'FaceColor', 'r');    
    set(h(2), 'FaceColor', 'b');
    ylabel('# Obs.');
%     legend({['Early ' strRange], ['Late ' strRange]}); xlabel('X'); ylabel('# Obs.');
    set(gca, 'XLim', xi([1 end]));
% end