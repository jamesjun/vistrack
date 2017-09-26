function plotKS2(data1, data2)
% two population comparison
%based on Loren Frank's version
confint = 95;
vrQt = 0:.001:1;
N = sqrt(numel(data1)^2 + numel(data2)^2);
vrQd = quantile(data1, vrQt); %data quantile

KSSorted = zeros(size(vrQd));
for i=1:numel(vrQd)
    KSSorted(i) = mean(data2 <= vrQd(i));
end

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
Y1 = KSSorted - vrQt; Y2 = zeros(size(vrQt)); ylim1 = [-10, 10] * KSdist/sqrt(N);


plot( vrQd , Y1, 'k', 'LineWidth', 2.0);  
hold on
h = plot(vrQt, Y2, 'LineWidth', 1); 
set(h, 'Color', [.4 .4 .4]);
h = plot(vrQd, Y2+KSdist/sqrt(N), '--', vrQd, Y2-KSdist/sqrt(N), '--' ); 
set(h(1), 'Color', [.4 .4 .4]);
set(h(2), 'Color', [.4 .4 .4]);
axis([vrQd(1) vrQd(end) ylim1(1) ylim1(2)]);
xlabel('Q1'); ylabel('Q2 - Q1');

probOver = sum(abs(KSSorted - vrQt) > KSdist/sqrt(N)) / N;
fprintf('KS test: %0.2f%% is outside the CI\n', probOver*100);
