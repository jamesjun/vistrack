function [mlMask, regionStr] = getImageMask(img0, rangeLim, strLoc)
% Get distance from landmark shapes

pixpercm = 1053.28/(sqrt(2)*100);

if nargin < 3
    strLoc = 'CENTRE';
end

switch upper(strLoc)
    case 'CENTRE'
        xy0 = [787, 606];
        mlMaskA = makeMask(xy0, rangeLim(1)*pixpercm*2, img0, 'CIRCLE');
        mlMaskB = makeMask(xy0, rangeLim(2)*pixpercm*2, img0, 'CIRCLE');  
        mlMask = mlMaskB & ~mlMaskA;
    case 'FOOD' %food, radius=1cm
        xy0 = [789, 681];
        mlMaskA = makeMask(xy0, (rangeLim(1)+1)*pixpercm*2, img0, 'CIRCLE');
        mlMaskB = makeMask(xy0, (rangeLim(2)+1)*pixpercm*2, img0, 'CIRCLE');  
        mlMask = mlMaskB & ~mlMaskA;
    case 'LM1' %landmark 1 (small square)
        xy0 = [966 418];
        mlMaskA = makeMask(xy0, 2.216*2.54*pixpercm, img0, 'SQUARE', rangeLim(1)*pixpercm);
        mlMaskB = makeMask(xy0, 2.216*2.54*pixpercm, img0, 'SQUARE', rangeLim(2)*pixpercm);  
        mlMask = mlMaskB & ~mlMaskA;
    case 'LM2' %landmark 2 (large square)
        xy0 = [975 790];
        mlMaskA = makeMask(xy0, 3.545*2.54*pixpercm, img0, 'SQUARE', rangeLim(1)*pixpercm);
        mlMaskB = makeMask(xy0, 3.545*2.54*pixpercm, img0, 'SQUARE', rangeLim(2)*pixpercm); 
        mlMask = mlMaskB & ~mlMaskA;
    case 'LM3' %landmark 3 (large circle)
        xy0 = [604 799];        
        mlMaskA = makeMask(xy0, (4*2.54+rangeLim(1)*2)* pixpercm, img0, 'CIRCLE');
        mlMaskB = makeMask(xy0, (4*2.54+rangeLim(2)*2)* pixpercm, img0, 'CIRCLE');
        mlMask = mlMaskB & ~mlMaskA;
    case 'LM4' %landmark 4 (small circle)
        xy0 = [600 428];
        mlMaskA = makeMask(xy0, (3*2.54+rangeLim(1)*2)* pixpercm, img0, 'CIRCLE');
        mlMaskB = makeMask(xy0, (3*2.54+rangeLim(2)*2)* pixpercm, img0, 'CIRCLE');
        mlMask = mlMaskB & ~mlMaskA;
    case 'LM*'
        mlMask = getImageMask(img0, rangeLim, 'LM1') | ...
                 getImageMask(img0, rangeLim, 'LM2') | ...
                 getImageMask(img0, rangeLim, 'LM3') | ...
                 getImageMask(img0, rangeLim, 'LM4');
    case 'LM*F'
        mlMask = makeMask([970 420], 2.216*2.54*pixpercm, img0, 'SQUARE', rangeLim(2)*pixpercm) | ...
                 makeMask([972 791], 3.545*2.54*pixpercm, img0, 'SQUARE', rangeLim(2)*pixpercm) | ...
                 makeMask([599 792], (4*2.54+rangeLim(2)*2)* pixpercm, img0, 'CIRCLE') | ...
                 makeMask([601 421], (3*2.54+rangeLim(2)*2)* pixpercm, img0, 'CIRCLE');             
    otherwise
        error(sprintf('%s not recognized.', strLoc));
end

% if ~strcmpi(strLoc, 'LM*')
%     [h, w] = size(img0);
%     [X, Y] = meshgrid(1:w, 1:h);
%     mlMask = false(size(img0));
%     vrR = sqrt((X - xy0(1)).^2 + (Y - xy0(2)).^2) / pixpercm;
%     mlMask(vrR >= rangeLim(1) & vrR < rangeLim(2)) = 1;
% end

regionStr = sprintf('%s %d~%d cm', strLoc, rangeLim(1), rangeLim(2));

if nargout == 0
    img0 = imadjust(img0);
    img0(~mlMask) = img0(~mlMask)/4;
    
    figure; imshow(img0);
    hold on; 
    if ~strcmpi(strLoc, 'LM*')
        plot(xy0(1), xy0(2), 'r.');
    end
end

end %func
