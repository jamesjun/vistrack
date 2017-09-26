function [tau, tau_sd] = calcCorrTau(vrX, strField)
% vrX: struct array

% if nargin < 2
    thresh = 1/exp(1);
% end

% array input
if iscell(vrX) || isstruct(vrX(1))
    vrTau = zeros(size(vrX));
    for i=1:numel(vrTau)
        if iscell(vrX)
            S = vrX{i};
        else
            S = vrX(i);
        end
        if isstruct(S)
            vrTau(i) = calcCorrTau(getfield(S, strField));
        else
            vrTau(i) = calcCorrTau(S);
        end
    end
    tau = mean(vrTau); %mean correlation coeff.    
    tau_sd = std(vrTau); %mean correlation coeff.    
    return;
end

vrC = xcorr(vrX - mean(vrX), 'coeff');
vrC = vrC(ceil(end/2):end);
tau = find(vrC < thresh, 1, 'first');
tau_sd = 0;
% figure; plot(vrC);

end