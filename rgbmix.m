%   mix the RGB to RGBbk in the masked area
function RGB = rgbmix(RGBbk, RGB, MASK, mode, mixRatio)
% RGB = rgbmix(RGBbk, RGB, MASK, 'mix', mixRatio)
% RGB = rgbmix(RGBbk, RGB, [], 'transparent', mixRatio)

if nargin<3, MASK = []; end
if nargin<4, mode = ''; end
if nargin<5, mixRatio = []; end

if isempty(mixRatio), mixRatio = .25; end
if isempty(mode)
    if ~isempty(MASK)
        mode = 'mix';
    else
        mode = 'transparent';
    end
end

if numel(size(RGBbk)) == 2 %gray scale
    if isa(RGBbk, 'double')
        RGBbk = uint8(RGBbk/max(RGBbk(:)));
    end
    RGBbk = imgray2rgb(RGBbk, [0 255], 'gray');
end
if numel(size(RGB)) == 2 %gray scale
    if ~isempty(MASK)
        RGB(~MASK) = 0;
    end
    if isa(RGB, 'double')
        RGB = uint8(RGB/max(RGB(:))*255);
    end
    RGB = imgray2rgb(RGB, [0 255], 'jet');
end

% R = RGBbk(:,:,1);
% G = RGBbk(:,:,1);
% B = RGBbk(:,:,1);

switch mode
    case 'mix'
        for iColor = 1:3
            mr_ = uint8(RGB(:,:,iColor)*mixRatio + RGBbk(:,:,iColor)*(1-mixRatio));            
            if isempty(MASK)
                RGB(:,:,iColor) = mr_;
            else
                mr1_ = RGB(:,:,iColor);
                mr1_(MASK) = mr_(MASK);
                RGB(:,:,iColor) = mr1_;
            end
        end
    case 'transparent'
        %mask is not used, instead RGBbk and RGB are simply added 50/50
        for iColor = 1:3
            RGB(:,:,iColor) = uint8(RGB(:,:,iColor)*mixRatio + RGBbk(:,:,iColor)*(1-mixRatio));
        end
    otherwise
        error('rgbmix invalid mode');
end
end %func