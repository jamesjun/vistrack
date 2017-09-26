function plotColorbar(dimm, vrRateSrt, vrQuantSrt)
% figure(hfig);

colormap jet;
rect = [dimm(2)-140, 50, 50, dimm(1)-100];
vrX = [rect(1); rect(1) + rect(3) - 1; rect(1) + rect(3) - 1; rect(1)];
vrY = [rect(2); rect(2); rect(2) + rect(4) - 1; rect(2) + rect(4) - 1];
patch(vrX, vrY, [256;256;1;1]);

nlabels = 4;
vnY = [(rect(2) + rect(4) - 1):-1:rect(2)];
% vrQuantile = [1:256]/256;
viTick = round((0:nlabels)/nlabels*256);
viTick = min(max(viTick,1),256);
vrQuantile1 = vrQuantSrt(viTick);
vrRateSrt1 = vrRateSrt(viTick);
ny = numel(vnY);
viY = round((0:nlabels)/nlabels * ny);
viY = min(max(viY,1),ny);
vnY1 = vnY(viY);
text(rect(1)-50, 550, 'EOD (Hz)',...
     'FontSize', 12, 'Color', [1 1 1], 'Rotation', 90);
text(rect(1)+85, 550, 'Quantile',...
     'FontSize', 12, 'Color', [1 1 1], 'Rotation', 90);

for i=1:(nlabels+1)
    text(rect(1)-73, vnY1(i), ...
         sprintf('%0.1f', vrRateSrt1(i)) ,...
         'FontSize', 12, 'Color', [1 1 1]);
     
    text(rect(1)+61, vnY1(i), ...
         sprintf('%0.2f', vrQuantile1(i)) ,...
         'FontSize', 12, 'Color', [1 1 1]);   
     
    line([rect(1)-2, rect(1) + rect(3)+2], vnY1(i)*[1 1],...
            'Color', [1 1 1]);
end

% hrect = rectangle('Position', [], 'FaceColor','r');