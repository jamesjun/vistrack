function RGB = gray2rgb(im, intlim)
if nargin < 2
    intlim = [0 255];
end
    im = imadjust(im, intlim/255, [0,1]);
    RGB = repmat(im, [1,1 , 3]);
end