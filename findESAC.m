function [vrESAC, viESAC, vrDur] = findESAC(EODA)

viUp = find(diff(EODA>0)>0)+1;
vrESAC = zeros(size(viUp));
viESAC = zeros(size(viUp));
vrDur = zeros(size(viUp));
for i=1:numel(viUp)    
    iStart = viUp(i);
    iEnd = find(EODA(iStart:end) <= 0, 1, 'first') - 1 + iStart;
    if ~isempty(iEnd)
        [vrESAC(i), idx] = max(EODA(iStart:iEnd));
        viESAC(i) = idx - 1 + iStart;
        vrDur(i) = iEnd - iStart;
    else
        vrESAC = vrESAC(1:i-1);
        viESAC = viESAC(1:i-1);
        return;
    end
end