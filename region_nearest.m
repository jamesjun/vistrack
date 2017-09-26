function [iRegion] = region_nearest(regions, xy0)
% get the region nearest
if numel(regions)==1, iRegion = 1; return; end
mrCentroid_region = cell2mat(arrayfun(@(x)x.Centroid, regions, 'UniformOutput', 0));
[~, iRegion] = min(pdist2(mrCentroid_region, xy0));
end