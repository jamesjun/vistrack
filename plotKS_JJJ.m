function plotKS_JJJ(data, PD)
%PD: probability distribution object
%based on Loren Frank's version
confint = 99;
vrQt = 0:.001:1;
N = numel(data);
vrQd = quantile(data, vrQt); %data quantile
KSSorted = cdf(PD, vrQd);

if (confint == 95)
    KSdist = 1.36;
elseif (confint == 99)
    KSdist = 1.63;
else
    error('Invalid confidence bound');
end

%obtain quantile from PD object
% ([1:N]-.5)/N

% Y1 = vrQT; Y2 = vrQT; ylim1 = [0 1];
Y1 = KSSorted - vrQt; Y2 = zeros(size(vrQt)); ylim1 = [-3, 3] * KSdist/sqrt(N);


plot( KSSorted , Y1, 'k', 'LineWidth', 2.0);  
hold on
h = plot(vrQt, Y2, 'LineWidth', 1); 
set(h, 'Color', [.4 .4 .4]);
h = plot(vrQt, Y2+KSdist/sqrt(N), '--', vrQt, Y2-KSdist/sqrt(N), '--' ); 
set(h(1), 'Color', [.4 .4 .4]);
set(h(2), 'Color', [.4 .4 .4]);
axis([0 1 ylim1(1) ylim1(2)]);
xlabel('Theoretical CDF'); ylabel('Empirical CDF');

probOver = sum(abs(KSSorted - vrQt) > KSdist/sqrt(N)) / N;
fprintf('KS test: %0.2f%% is outside the CI\n', probOver*100);
