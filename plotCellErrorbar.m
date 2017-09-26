function mrY = plotCellErrorbar(cellMat, strVar, fun1, kcorr)
%n by 3 (val, low, high)
%fun1: override function

if nargin == 1
    mrY = plotCellErrorbar(cellMat, '', [],1);
    return;
end

[m, n] = size(cellMat);
nBoot = 250;

bootMean = @(x)[mean(x); bootci(nBoot, {@(y)mean(y), x}, 'type', 'cper')];

vrX = 1:m;
mrY = zeros(m, n);
mrH = zeros(m, n);
mrL = zeros(m, n);

if nargin < 3
    fun1 = @(x)x;
elseif isempty(fun1)
    fun1 = @(x)x;
end
if nargin < 4
    fcorr = 1; %calculate on-fly. be careful of multiple trial types
else
    fcorr = 0; %use given kcorr value
end

for i=1:m
    for j=1:n
        vr = cellMat{i, j};        
        if isstruct(vr)
            vr = getfield(vr, strVar);
        end
%         vr = fun1(vr(~isnan(vr) & ~isinf(vr)));

        vec = sort(fun1(bootMean(vr)));
        if fcorr
            kcorr = calcCorrTau(vr); %correlation correction  
        end

        mrY(i,j) = vec(2);
        mrL(i,j) = (vec(2) - vec(1)) * sqrt(kcorr);
        mrH(i,j) = (vec(3) - vec(2)) * sqrt(kcorr);
    end
end
if n>1
    vrX1 = linspace(-.1*n, .1*n, n) * .75; %((1:n) - 1 - round(n/2)) / n * .65;
else
    vrX1 = 0;
end
h = bar(vrX, mrY, .75); hold on;

for i=1:m        
    errorbar(i + vrX1, mrY(i,:), mrL(i,:), mrH(i,:), 'r.'); hold on;
end
set(gca, 'XLim', [.5 m+.5]);
set(gca, 'XTick', vrX);

%color
if numel(h) == 3
    set(h(1), 'FaceColor', 'r');
    set(h(2), 'FaceColor', 'b');
    set(h(3), 'FaceColor', 'g');
elseif numel(h) == 2
    set(h(1), 'FaceColor', 'k');
    set(h(2), 'FaceColor', [.5 .5 .5]);
end

