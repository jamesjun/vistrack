function [pValue, KSstatistic] = pval_kstest2(x1, x2, n1, n2)

tail = 0; %two tailed
%
% Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
%

binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

%
% Compute the test statistic of interest.
%

switch tail
   case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
      deltaCDF  =  abs(sampleCDF1 - sampleCDF2);

   case -1      %  1-sided test: T = max[F2(x) - F1(x)].
      deltaCDF  =  sampleCDF2 - sampleCDF1;

   case  1      %  1-sided test: T = max[F1(x) - F2(x)].
      deltaCDF  =  sampleCDF1 - sampleCDF2;
end

KSstatistic   =  max(deltaCDF);

%
% Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.
%

n      =  n1 * n2 /(n1 + n2);
lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);

if tail ~= 0        % 1-sided test.

   pValue  =  exp(-2 * lambda * lambda);

else                % 2-sided test (default).
%
%  Use the asymptotic Q-function to approximate the 2-sided P-value.
%
   j       =  (1:101)';
   pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
   pValue  =  min(max(pValue, 0), 1);

end
