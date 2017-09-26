function [vrYf, slope, intercept] = fitPoly1(vrX, vrY)

if min(size(vrX)) == 1
    vrX = vrX(:);
    vrY = vrY(:);
end
n = size(vrX,1);

mrX = [ones(n,1), vrX];
M = inv(mrX' * mrX) * mrX';

B = M * vrY;
slope = B(2:end);
intercept = B(1);
vrYf = vrX * slope + intercept;

if nargout == 0
    figure; plot(vrX, vrY, 'b', vrX, vrYf, 'r'); 
    legend('raw', 'fitted');
end