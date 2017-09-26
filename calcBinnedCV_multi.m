function [X_cv] = calcBinnedCV_multi( X, nbin )
%CALCBINNEDMEAN Summary of this function goes here
%   X: nch by nbins
%   length(Y) = floor(numel(X)/binsize)
%   length(Y1) = length(X), Y1 repeats and preserves the length, returns a
%   column vector for Y and Y1 for any vector format of X

X_cv = calcBinnedCV(X(:,1), nbin);
for i=2:size(X,2)
    X_cv = X_cv + calcBinnedCV(X(:,i), nbin);
end