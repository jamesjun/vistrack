function [mrColor, vrRateSrt, vrQuantSrt] = quantile2color(RC)
RC=RC(:);
nIndex = 256; %256 level color
pwr = 3;
mrJet = jet(nIndex);
n = numel(RC);

% calculate rank
viRank = zeros(1, n);
for i=1:numel(RC)
    viRank(i) = sum(RC < RC(i));
end

% scale color
% range from -1 to 1 and apply gamma power
viColor = (viRank - n/2) / (n/2);
viColor = sign(viColor) .* abs(viColor) .^ pwr;
viColor = (viColor + 1) / 2 * nIndex; % range from 1~nIndex
viColor1 = round(viColor);
viColor1 = max(min(viColor1, nIndex), 1);
mrColor = mrJet(viColor1,:);

%colorbar scale
[viColorSrt, IX] = sort(viColor, 'ascend');
RCSrt = RC(IX);
vrQuantSrt = viRank(IX)/nIndex;
vrRateSrt = interp1(viColorSrt, RCSrt, 1:nIndex, 'spline', 'extrap');

end