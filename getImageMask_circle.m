function [mlMask, regionStr] = getImageMask_circle(img0, rangeLim, strLoc)
% Get circular mask based on distance only

pixpercm = 1053.28/(sqrt(2)*100);

if nargin < 3
    strLoc = 'CENTRE';
end

switch upper(strLoc)
    case 'CENTRE'
        xy0 = [787 606];
    case 'FOOD'
        xy0 = [789  681]; %food
    case 'LM1'
        xy0 = [966 418]; %landmark 1 (small square)
    case 'LM2'
        xy0 = [975 790]; %landmark 1 (small square)
    case 'LM3'
        xy0 = [604 799]; %landmark 1 (small square)
    case 'LM4'
        xy0 = [600 428]; %landmark 1 (small square)
    case 'LM*'
        mlMask = getImageMask(img0, rangeLim, 'LM1') | ...
                 getImageMask(img0, rangeLim, 'LM2') | ...
                 getImageMask(img0, rangeLim, 'LM3') | ...
                 getImageMask(img0, rangeLim, 'LM4');
    otherwise
        error(sprintf('%s not recognized.', strLoc));
end

if ~strcmpi(strLoc, 'LM*')
    [h, w] = size(img0);
    [X, Y] = meshgrid(1:w, 1:h);
    mlMask = false(size(img0));
    vrR = sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;
    mlMask(vrR >= rangeLim(1) & vrR < rangeLim(2)) = 1;
end

regionStr = sprintf('%s %d~%d cm', strLoc, rangeLim(1), rangeLim(2));

if nargout == 0
    img0 = imadjust(img0);
    img0(~mlMask) = img0(~mlMask)/2;
    
    figure; imshow(img0);
    hold on; 
    if ~strcmpi(strLoc, 'LM*')
        plot(xy0(1), xy0(2), 'r.');
    end
end

end