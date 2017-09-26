function RGB = imgray2rgb(I, inputrange, cmap)
%   imgray2rgb  converts image to RGB scaled image, unit8
%   JJJ function

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

if ~exist('cmap');
    cmap = 'Jet';
end

if strcmp(class(I), 'uint8')
    if ~isempty(inputrange)        
        I = imadjust(I, double(inputrange)/255, [0 1]);
    end
    I = uint8(I/4); %index ranges only from 1 to 64 whereas uint8 ranges 0 5o 255
    RGB = ind2rgb8(I, colormap(cmap));
else
    I = double(I);
    theCmap = colormap(jet(255));
    I = uint8((I-inputrange(1)) / (inputrange(2)-inputrange(1)) * size(theCmap, 1));
    RGB = ind2rgb8(I, theCmap);
end