figure; 
colormap jet;
rect = [500, 1, 10, 300];
vrX = [rect(1); rect(1) + rect(3) - 1; rect(1) + rect(3) - 1; rect(1)];
vrY = [rect(2); rect(2); rect(2) + rect(4) - 1; rect(2) + rect(4) - 1];
patch(vrX, vrY, [1;1;256;256]);


% patch([.25;.75;.75;.25], [0;0;.5;.5], [1 1 1]);
% axis([-2, 2, -2, 2]);

%      clear cdata 
% set(gca,'CLim',[0 40])
% cdata = [15 30 25 2 ...
% 60 12 23 40 13 26 24]';
% set(p,'FaceColor','interp',...
% 'FaceVertexCData',cdata,...
% 'EdgeColor','flat',...
% 'LineWidth',0,...
% 'CDataMapping','scaled')