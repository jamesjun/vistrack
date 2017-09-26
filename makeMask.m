function mlMask = makeMask(xy0, d1, img0, strShape, r)
% d1: diameter
% r: range expansion

if nargin < 5
    r = 0;
end

d1 = round(d1);
if d1 < 1
    mlMask = false(size(img0)); %none included
    return;
end
if nargin < 4
    strShape = 'CIRCLE';
end

fig = figure; 
warning off; 
image(false(size(img0)));
switch upper(strShape)
    case 'CIRCLE'
    h = imellipse(gca, [xy0(1)-d1/2, xy0(2)-d1/2, d1, d1]); %[x y w h]
    case 'SQUARE'
    h = imrect(gca, [xy0(1)-d1/2, xy0(2)-d1/2, d1, d1]); %[x y w h]
    case 'RECT'
    h = imrect(gca, [xy0(1)-d1(1)/2, xy0(2)-d1(2)/2, d1(1), d1(2)]); %[x y w h]    
end
mlMask = createMask(h);

% make a round square
if r >= 1 && strcmpi(strShape, 'SQUARE')
    mlMask1 = makeMask(xy0 + [+1,+1]*d1/2, 2*r, img0, 'CIRCLE');
    mlMask2 = makeMask(xy0 + [+1,-1]*d1/2, 2*r, img0, 'CIRCLE');
    mlMask3 = makeMask(xy0 + [-1,+1]*d1/2, 2*r, img0, 'CIRCLE');
    mlMask4 = makeMask(xy0 + [-1,-1]*d1/2, 2*r, img0, 'CIRCLE');
    mlMask5 = makeMask(xy0, [d1, d1+2*r], img0, 'RECT');
    mlMask6 = makeMask(xy0, [d1+2*r, d1], img0, 'RECT');
    mlMask = mlMask1 | mlMask2 | mlMask3 | mlMask4 | mlMask5 | mlMask6;
end
    
%expand

if nargout > 0
    close(fig);
else
    img1 = imadjust(img0);
    img1(mlMask) = 0;
    imshow(img1); 
end
end