function [BW blobsize nAreas] = bwgetlargestblob(BW, blobsizelim)

if nargin == 1
    blobsizelim = [-1 inf];
end
nAreas = 0;
[L, nL] = bwlabel(BW, 4);
if nL == 1    
    blobsize = sum(BW(:));
    if (blobsize > blobsizelim(2)) || (blobsize < blobsizelim(1)) %out of range
        BW = false(size(BW));
        blobsize = 0;
    else
        nAreas = 1;
    end
    return;
end

%search for the largest blob
blobsize=0;
imax = nan;
for i=1:nL
    area = sum(sum(L==i));
    if area < blobsizelim(2) && area > blobsizelim(1) %out of range
        nAreas = nAreas + 1;
        if area > blobsize
            blobsize = area;
            imax = i;
        end
    end
end
if ~isnan(imax)
    BW = (L==imax);
else
    BW = false(size(BW));
end