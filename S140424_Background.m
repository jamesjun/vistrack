% show image and fish approaching and moving away

load img0_P.mat;
rectCrop = [493 1083 312 902];

img0 = imadjust(img0_P);

img0 = imrotate(img0, -1.1590, 'nearest', 'crop');
figure; imshow(img0);
axis(rectCrop);

%%
fname = 'C:\expr\Raw data\E13A4p.wmv';
vid = VideoReader(fname);
mov = read(vid, [1 1500]);
implay(mov);

%%
img1 = mov(:,:,1,1180);
img2 = mov(:,:,1,1200);

img1 = imadjust(img1);
img2 = imadjust(img2);

img = min(img1, img2);
figure; imshow(img);


%%
fname = 'C:\expr\Raw data\E13A4p.wmv';
vid = VideoReader(fname);
mov = read(vid, [1 1500]);
implay(mov);


%%
iFrame = 1231;
img0 = mov(:,:,1,iFrame);
img0 = imadjust(img0);

img0 = imrotate(img0, -1.1590, 'nearest', 'crop');
figure; imshow(img0);
axis(rectCrop);