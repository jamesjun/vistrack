function RGB = imgray2rgb(I, inputrange, vcColorMap)
%   imgray2rgb  converts image to RGB scaled image, unit8
%   JJJ function
if nargin<3, vcColorMap = 'jet'; end  % PARULA, HSV, HOT, PINK

if ~exist('inputrange')
    if strcmp(class(I), 'uint8')
        inputrange = [0 255];
    else
        inputrange = [min(I(:)) max(I(:))];
    end
else
    if isempty(inputrange)
        if ~strcmp(class(I), 'uint8')
            inputrange = [min(I(:)) max(I(:))]; 
        end
    end
end

if strcmp(class(I), 'uint8')
    if ~isempty(inputrange)        
        I = imadjust(I, double(inputrange)/255, [0 1]);
    end
    I = uint8(I);
    RGB = ind2rgb_(I, vcColorMap);
else
    I = double(I);
    I = uint8((I-inputrange(1)) / (inputrange(2)-inputrange(1)) * size(theCmap, 1));
    RGB = ind2rgb_(I, vcColorMap);
end
end %func


%--------------------------------------------------------------------------
function [RGB, mrCmap] = ind2rgb_(miImg, vcColorMap)
% I: int8, vcColorMap: char
eval(sprintf('mrCmap = uint8(%s(256)*255);', lower(vcColorMap)));
RGB = zeros([size(miImg), 3], 'uint8');
miImg = miImg + 1; % 0 base to 1 base
for iColor = 1:3
    viMap_ = mrCmap(:,iColor);
    RGB(:,:,iColor) = viMap_(miImg);
end
end %func