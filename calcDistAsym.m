function delMed = calcDistAsym(vrY)

% delMed = skewness(vrY); return;

sd = std(vrY);
delMed = sum(vrY>sd) / sum(vrY<-sd) - 1;

%EODA+ median + EODA- median (negative)
% 
% vrL = vrY>=0;
% delMed = std(vrY(vrL)) - std(vrY(~vrL));
% 

if nargout==0
    [~,xi,bw] = ksdensity(vrY);
    vrX = xi(1):bw:xi(end);
    figure;
    ksdensity(vrY, vrX, 'bandwidth', bw); hold on;
    plot(sd*[1 1], get(gca, 'YLim'), 'r-');
    plot(sd*[-1 -1], get(gca, 'YLim'), 'r-');
end
