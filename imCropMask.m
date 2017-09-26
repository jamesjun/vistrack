function mrImgCrop = imCropMask(mrImg, mlMask)

rp = regionprops(mlMask, 'BoundingBox');
bb = round(rp.BoundingBox);
mrImgCrop = mrImg(bb(2):bb(2)+bb(4), bb(1):bb(1)+bb(3));
