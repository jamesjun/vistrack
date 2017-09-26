function vrN = normIPI(vrI)
%normalize IPI using caputi's method
% 1-I/I0

% Mode normalization
[f, xi, bw] = ksdensity(vrI);
[~,imax] = max(f);
I0 = xi(imax); 

% I0 = mean(vrI);
% I0 = quantile(vrI, .5);

% vrN = zscore(vrI);
vrN = vrI/I0;

if nargout == 0
    figure;
    subplot 211;
    ksdensity(vrI); hold on;
%     cdfplot(vrI);
    plot(I0*[1 1], get(gca, 'YLim'));

    subplot 212;
    plot(vrN); hold on;
    plot(get(gca, 'XLim'), [1 1]);
end