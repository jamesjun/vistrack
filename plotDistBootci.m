function [xi, vrF, vrF_L, vrF_H] = plotDistBootci(vrX)

[vrF,xi,bw] = ksdensity(vrX, 'function', 'survivor'); 
mrFci = bootci(100, {@(x)ksdensity(x, xi, 'function','survivor','bandwidth',bw),vrX}, 'type', 'cper');

vrF_L = mrFci(1,:);
vrF_H = mrFci(2,:);

vrF = vrF(:);
vrF_L = vrF_L(:);
vrF_H = vrF_H(:);

if nargout == 0
figure; hold on;
plot(xi, vrF);
plot(xi, vrF_low);
plot(xi, mrFci(2,:));
end