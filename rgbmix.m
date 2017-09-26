function RGB = rgbmix(RGBbk, RGB, MASK, mode)
%   mix the RGB to RGBbk in the masked area
mixRatio = .25;

if nargin <3
    mode = 'transparent';
elseif nargin < 4
    mode = 'mix';
end

if numel(size(RGBbk)) == 2 %gray scale
    if isa(RGBbk, 'double')
        RGBbk = uint8(RGBbk/max(RGBbk(:)));
    end
    RGBbk = imgray2rgb(RGBbk, [0 255], 'Gray');
end
if numel(size(RGB)) == 2 %gray scale
    if nargin >= 3
        RGB(~MASK) = 0;
    end
    if isa(RGB, 'double')
        RGB = uint8(RGB/max(RGB(:))*255);
    end
    RGB = imgray2rgb(RGB, [0 255], 'Jet');
end

R = RGBbk(:,:,1);
G = RGBbk(:,:,1);
B = RGBbk(:,:,1);

switch mode
    case 'mix'
        Rmix = uint8((RGB(:,:,1)*mixRatio + R*(1-mixRatio)));
        Gmix = uint8((RGB(:,:,2)*mixRatio + G*(1-mixRatio)));
        Bmix = uint8((RGB(:,:,3)*mixRatio + B*(1-mixRatio)));

        R(MASK) = Rmix(MASK);
        G(MASK) = Gmix(MASK);
        B(MASK) = Bmix(MASK);
    case 'transparent'
        %mask is not used, instead RGBbk and RGB are simply added 50/50
        R = uint8((RGB(:,:,1)*mixRatio + R*(1-mixRatio)));
        G = uint8((RGB(:,:,2)*mixRatio + G*(1-mixRatio)));
        B = uint8((RGB(:,:,3)*mixRatio + B*(1-mixRatio)));          
    otherwise
        error('rgbmix invalid mode');
end

RGB(:,:,1) = R;
RGB(:,:,2) = G;
RGB(:,:,3) = B;