function [Corners, MASK] = getBoundingBoxPos(CenterPos, imagesize, winsize)
%GETBOUNDINGBOXPOS Calculate the coordinates of the cropping window
%[Input]
% CenterPos: [x, y]. Center of the cropping window
% imagesize: [height, width]. Size of the image
% winsize: [height, width]. Size of the cropping window
%[Output]
% Corners: [x1,x2,y1,y2]. top left and bottom right corners
% MASK: binary image mask having the same size as the image.

% Default parameters
if nargin < 3
    winsize = [192 192];
end    
if nargin < 2
    imgheight = 480;
    imgwidth = 640;
else
    imgheight = imagesize(1);
    imgwidth = imagesize(2);
end

% format inputs
winheight = winsize(1);
winwidth = winsize(2);

% Calculate corner coordinates
x1 = round(CenterPos(1) - winwidth/2);
y1 = round(CenterPos(2)- winheight/2);
if x1 < 1
    x1 = 1;
end
x2 = x1 + winwidth - 1;
if x2 > imgwidth
    x2 = imgwidth;
    x1 = x2 - winwidth + 1;
end
if y1 < 1
    y1 = 1;
end
y2 = y1 + winheight - 1;
if y2 > imgheight
    y2 = imgheight;
    y1 = y2 - winheight + 1;
end
Corners = [x1, x2, y1, y2];

% Create an image mask
if nargout == 2
    MASK = false(imgheight, imgwidth);
    MASK(y1:y2, x1:x2) = 1;
end
end %func