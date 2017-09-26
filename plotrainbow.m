function mrC = plotrainbow(X, Y, markerSize)
% time-plot rainbow
if nargin < 3
    markerSize = 4;
end

mrC = colormap(cool(numel(X)));
for i=1:numel(X)
    plot(X(i), Y(i), '.', 'color', mrC(i,:)); %, 'MarkerSize', markerSize);
end

end