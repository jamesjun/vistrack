function plotBarMuSem(csZ, XTickLabel, csColor)

mrMu = zeros(size(csZ));
mrSem = zeros(size(csZ));
for i=1:size(csZ,1)
    for j=1:size(csZ,2)
        vrZ = csZ{i,j};
        mrMu(i,j) = mean(vrZ(:));
%         mrSem(i,j) = std(vrZ(:))/numel(vrZ);
        mrSem(i,j) = std(vrZ(:));
    end
end

plotBarError(mrMu, mrSem, csColor, XTickLabel);

