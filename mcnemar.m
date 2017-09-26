function pval = mcnemar(vl1, vl2)
N10 = sum(logical(vl1) & ~logical(vl2));
N01 = sum(~logical(vl1) & logical(vl2));

mcnemar = (abs(N10-N01) - 1)^2/(N10+N01);
pval = 1 - chi2cdf(mcnemar,1);