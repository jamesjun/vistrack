function mrRGB = maskRGB(mrRGB, MASK, val)
if nargin < 3
    val = 0;
end

mrR = mrRGB(:,:,1);
mrG = mrRGB(:,:,2);
mrB = mrRGB(:,:,3);

mrR(MASK) = mrR(MASK) * val;
mrG(MASK) = mrG(MASK) * val;
mrB(MASK) = mrB(MASK) * val;

mrRGB(:,:,1) = mrR;
mrRGB(:,:,2) = mrG;
mrRGB(:,:,3) = mrB;