%% image analysis code

figure; imshow(img2)
img = read(aviobj, FLIM(1)+100); figure; imshow(img);
img = read(aviobj, FLIM(1)+1000); figure; imshow(img);
img = read(aviobj, FLIM(1)+5000); figure; imshow(img);
img = read(aviobj, FLIM(1)+10000); figure; imshow(img);
img = read(aviobj, FLIM(1)+3000); figure; imshow(img);
img = read(aviobj, FLIM(1)+4000); figure; imshow(img);
figure; imagesc(img0*.9 - img)
size(img0)
size(img)
img = img(:,:,1);
figure; imshow(img0*.9 - img)
figure; imshow(img0*.9 - img,[0 25])
figure; imshow(img);
%-- 2/13/2013 11:35 AM --%