function RGB = ind2rgb8(X, CMAP)
%IND2RGB8 Convert indexed image to uint8 RGB image
%
%   RGB = IND2RGB8(X,CMAP) creates an RGB image of class uint8.  X must be
%   uint8, uint16, or double, and CMAP must be a valid MATLAB colormap.
%
%   Example 
%   -------
%      % Convert the 'concord_ortho_e.tif' image to RGB.
%      [X, cmap] = imread('concord_ortho_e.tif');
%      RGB = ind2rgb8(X, cmap);
%      R = worldfileread('concord_ortho_e.tfw');
%      mapshow(RGB, R);
%
%   See also IND2RGB.

% Copyright 1996-2007 The MathWorks, Inc.

%RGB = ind2rgb8c(X, CMAP);
try
    RGB = zeros([size(X,1)*size(X,2), 3], 'uint8');
    CMAP = uint8(CMAP * 255);
    nColors = size(CMAP, 1);
    X = min(max(X, 1), nColors);
    RGB(:,1) = CMAP(X,1);
    RGB(:,2) = CMAP(X,2);
    RGB(:,3) = CMAP(X,3);
    RGB = reshape(RGB, [size(X,1), size(X,2), 3]);
catch
	disp('err')   
end