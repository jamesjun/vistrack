function [X_cv X_mean X_std] = calcBinnedCV( X, binsize )
%CALCBINNEDMEAN Summary of this function goes here
%   length(Y) = floor(numel(X)/binsize)
%   length(Y1) = length(X), Y1 repeats and preserves the length, returns a
%   column vector for Y and Y1 for any vector format of X

nbins = floor(numel(X)/binsize);

X1 = reshape(X(1:binsize*nbins),binsize ,[]);

X_mean = mean(X1,1);
X_std = std(X1,[],1);
X_cv = X_std ./ X_mean; 


end