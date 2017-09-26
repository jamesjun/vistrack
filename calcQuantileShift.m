function Q = calcQuantileShift(X)
% gives high value for higher quantile shift

bin = 0:.001:1;
cdf = quantile(X, bin);

%remove repeated values
IKILL = [];
for i=2:numel(cdf)
    if (cdf(i-1) == cdf(i))
        IKILL(end+1) = i;
    end
end
cdf(IKILL) = [];
bin(IKILL) = [];

Q = interp1(cdf, bin, X, 'linear');
Q = Q-.5;
% figure;plot(cdf,bin);

% nwin=10;
% Q = filtfilt(ones(nwin,1)/nwin, 1, Q);