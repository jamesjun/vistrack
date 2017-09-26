function [img0, MASK, xyc, vec, xy0] = makeBackground(img1, img2)
%xyc: center of the fish's initial position
%xy0: center of the aquarium
%vec: initial heading direction (unit vec)
% persistent MaskPos xy0_init;
% default locations
backgroundMode = 3; 
%1: img2 as background
%2: img1 as major, img2 as minor
%3: img1 as minor, img2 as major

try
    load makeBackground_cache;
catch
    MaskPos = [205, 13, 1166, 1163];
    xy0_init = [794, 599];
end

hfig = figure;
set(hfig, 'OuterPosition', get(0, 'ScreenSize')); drawnow;
subplot 121; imshow(img1);  title('First frame');
subplot 122; imshow(img2);  title('Second frame');

% Crop out a fish from the first image
subplot 121;
uiwait(msgbox('Drag a box around the fish on the LEFT image and double-click', 'modal'));
title('Draw a box around and double-click');
h  = imrect;
position = wait(h);
MASK = createMask(h);
xc=round(position(1)+position(3)/2);
yc=round(position(2)+position(4)/2);
xyc = [xc yc];

% Define the head orientation
% uiwait(msgbox('Click on the head', 'modal'));
title('Click on the head');
xyh = ginput(1);
vec = xyh - xyc;
vec = vec ./ norm(vec);

switch backgroundMode
    case 1 %img2 as background
    img0 = img2;
    
    case 2 %img1 as major, img2 as minor
%     img0 = img1;    
%     img0(MASK) = img2(MASK);   
    img0 = mask_copy(img1, img2, MASK);
    
    case 3 %img1 as minor, img2 as major
    subplot 122;
    uiwait(msgbox('Drag a box around the fish on the RIGHT image and double-click', 'modal'));
    title('Draw a box around and double-click');
    h  = imrect;
    position = wait(h);
    MASK = createMask(h);        
%     img0 = img2;
%     img0(MASK) = img1(MASK);     
    img0 = mask_copy(img2, img1, MASK);
end

% Draw a mask
figure(hfig); clf; imshow(img0);
uiwait(msgbox('Create a circular mask and double-click', 'modal'));
title('Create a circular mask and double-click');
if exist('MaskPos', 'var')
    h = imellipse(gca, MaskPos);
else
    h  = imellipse;
end
wait(h);
% getPosition(h)
MaskPos = getPosition(h);
MASK = createMask(h);
img0 = mask_copy(img0, [], ~MASK);
% img0(~MASK) = 0;

% Define the center of the aquarium
% uiwait(msgbox('Click at the center of the aquarium', 'modal'));
% xy0 = round(pos_ellipse([1 2]) + pos_ellipse([3 4])/2);
xy0 = xy0_init;
title('Click on the center of the aquarium');
fAskUser = 0;
if exist('xy0', 'var')
    figure(hfig); imshow(img0); %display
    hold on; plot(xy0(1), xy0(2), 'r.');
    button = questdlg('Is the red dot centered?','Confirmation','Yes','No','Yes');
    if strcmp(button, 'Yes')
        try close(hfig), catch, end;
        return;
    else
        fAskUser = 1;
    end
else
    fAskUser = 1;
end

if fAskUser
    CenterPos = round([size(img1,2), size(img1,1)]/2);
    h = impoint(gca, CenterPos(1), CenterPos(2));
    xy0 = wait(h);
    xy0_init = xy0;
    
    % Final confirmation
    figure(hfig); imshow(img0); %display
    hold on; plot(xy0(1), xy0(2), 'r.');
    title('Composite background image');
    uiwait(msgbox('Composite background is created. Press OK to close.'));
    close(hfig);
end

save makeBackground_cache MaskPos xy0_init;
