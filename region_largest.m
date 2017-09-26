function [iRegion] = region_largest(regions)
% get the region nearest
if numel(regions)==1, iRegion = 1; return; end
[~, iRegion] = max([regions.Area]);
end