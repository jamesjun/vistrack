function mnVisit = calcVisitTime(vsTrialPool, img0)

vrX = poolVecFromStruct(vsTrialPool, 'vrX');
vrY = poolVecFromStruct(vsTrialPool, 'vrY');

H = fspecial('gaussian', [100 100], 20); %figure; imagesc(H);

mnVisit = zeros(size(img0));
mnVisit(sub2ind(size(img0), round(vrY), round(vrX))) = 1;
mnVisit = imfilter(mnVisit, H);


% if nargin < 3
%     nGrid = 20;
% end
% 
% vrX = poolVecFromStruct(vsTrialPool, 'vrX');
% vrY = poolVecFromStruct(vsTrialPool, 'vrY');
% 
% viX = ceil(vrX/nGrid);
% viY = ceil(vrY/nGrid);
% [h, w] = size(img0);
% h = h / nGrid;
% w = w / nGrid;
% 
% mnVisit = zeros(h, w);
% for iy=1:h
%     vlY = (viY == iy);
%     for ix=1:w
%         mnVisit(iy,ix) = sum(vlY & (viX == ix));
%     end
% end
% 
% if nargout == 0
%     mrVisit = uint8(mnVisit / max(mnVisit(:)) * 255);
%     
%     figure;
% %     imagesc(mrVisit);
%     imshow(rgbmix(img0, imresize(mrVisit, nGrid, 'nearest')));
% end
