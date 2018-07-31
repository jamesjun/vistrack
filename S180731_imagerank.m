% compute image difference
img_med = single(median(trImg0,3));
dimm0 = size(trImg0);
mrImg0 = single(reshape(trImg0, [], dimm0(3)));
mrImg_d = bsxfun(@minus, mrImg0, img_med(:)); % + bsxfun(@minus, img_med(:), mrImg0);
img_med2 = median(abs(mrImg_d),2);
figure; imshow(reshape(uint8(abs(mrImg_d(:,1))), dimm0(1:2)));


% color matching difference
% figure; plot(img_med(:), mrImg0(:,1), '.');
miImg0 = uint8(sortnorm(mrImg0')' * 255);
% [~,miImg0] = sort(mrImg0,2);
% miImg0(miImg0) = repmat((1:dimm0(3)), [dimm0(1)*dimm0(2), 1]);
% miImg0 = uint8(miImg0 / dimm0(3) * 255);
figure; imshow(reshape(miImg0(:,1), dimm0(1:2)));

img1 = mrImg0(:,1);
img0 = img_med(:);
imgd = reshape(uint8(abs(img1-img0)), dimm0(1:2));
imgr = reshape(img1 ./ img0, dimm0(1:2));
figure; imshow(imgr);
img2 = reshape(uint8(abs(img1 / nanmedian(imgr(:)) - img0)), dimm0(1:2));
figure; imshow(imadjust(img2));
figure; imshow(imadjust(imgd));
figure; imshow(reshape(uint8(img0), dimm0(1:2)));

a = single(mrImg0(:,1)) - single(img_med(:));
figure; imshow(imadjust(uint8(abs(reshape(a, dimm0(1:2))))));
ml = a~=0;
median(a(
int_med = median(a(ml));
a(ml) = a(ml) - int_med;
figure; imshow(imadjust(uint8(abs(reshape(a, dimm0(1:2))))));