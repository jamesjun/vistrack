function plotError(vrX, vrY, vrL, vrH, color)
% create 4x6 figure tile for direct export to illustrator

vrX = vrX(:);
vrY = vrY(:);
vrL = vrY - vrL(:);
vrH = vrY + vrH(:);

hold on;
% plot([vrX; flipud(vrX); vrX(1)], [vrL; flipud(vrH); vrL(1)], 'color', color);
fill([vrX; flipud(vrX)], [vrL; flipud(vrH)], color, 'EdgeColor', 'none'); %, 'FaceAlpha', .25);
% facealpha set will export as bitmap in eps
plot(vrX, vrY, 'color', color);