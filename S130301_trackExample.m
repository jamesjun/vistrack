
FLIM0 = [129 8230]; 
TLIM0 = [10 550]; 

vidFname = 'D:\\Expr\\121212_spontvid8\\121212_spontvid8ch000-1.wmv';
outFname = 'D:\\Expr\\121212_spontvid8\\VISTRACK_121212_000.mat';  

winlen = 128*2.5;
FPS_cam = 15.001551480959098;
SE = strel('disk',3,0); %use disk size 1 for 640 480
ThreshLim = [0 15];
IM_THRESH = 5;
IM_CONTRAST = .85;


disp('Loading the aviobj...'); 
tic; aviobj = VideoReader(vidFname); toc
disp('aviobj loaded')

%% pick first and last frame
if ~exist('FLIM0')
    implay(read(aviobj,[1 300]));
    uiwait(msgbox('Find first blink and close the movie'));
    ans = inputdlg('Frame Number', 'First frame',1,{'150'});
    FLIM0(1) = str2num(ans{1});

    FLIM0(2) = round(FLIM0(1) + FPS_cam * diff(TLIM0));
    FLIM1 = FLIM0(2) + [-50 50]; 
    FLIM1(2) = min(FLIM1(2), aviobj.NumberOfFrames);
    implay(read(aviobj,FLIM1));
    uiwait(msgbox('Find first blink and close the movie'));
    ans = inputdlg('Frame Number', 'Last frame',1,{'50'});
    temp = str2num(ans{1});
    FLIM0(2) = FLIM0(2)-50+temp-1;
end
if ~exist('FLIM')
    FLIM = FLIM0; %process all frames within the range
end
% determine time stamp
TC = interp1(FLIM0, TLIM0, FLIM(1):FLIM(2), 'linear');
TLIM = interp1(FLIM0, TLIM0, FLIM([1 end]), 'linear');
FPS = diff(FLIM0)/diff(TLIM0);
fprintf('FPS = %0.6f Hz, TLIM0=[%d, %d], FLIM0=[%d, %d]\n', FPS, TLIM0(1), TLIM0(2), FLIM0(1), FLIM0(2));


%% obtain and check background
img1 = read(aviobj, FLIM0(1)); img1=img1(:,:,1);
img2 = read(aviobj, FLIM0(2)); img2=img2(:,:,1);
% img2 = read(aviobj, 1000); img2=img2(:,:,1);

[img00, MASK, xy_init, vec0, xy0] = make_BG_JJJ(img1, img2);

%% obtain stretched background
intrng = stretchlim(img1, [0 .85]);
img1s = imadjust(img1, intrng);
img2s = imadjust(img2, intrng);
[img00, MASK, xy_init, vec0, xy0] = make_BG_JJJ(img1s, img2s);

%% subtract image
% implay(read(aviobj,[1 1000]));

iframe = 764;
xy0 = [805 1068];

img = read(aviobj, iframe); img=img(:,:,1);
imgs = imadjust(img, intrng);

% figure; imshow(img00)
figure; imshow(imgs);

xy1 = xy0 - 128;
xy2 = xy0 + +128 + 1;
imgsc = imgs(xy1(2):xy2(2), xy1(1):xy2(1)); figure; imshow(imgsc);
img00c = img00(xy1(2):xy2(2), xy1(1):xy2(1));   figure; imshow(img00c);
imgd = img00c*.85 - imgsc;  figure; imshow(imgd);
bw = imgd > 20;  figure; imshow(bw);






% bw1 = imclose(bw, ones(5, 5)); 
% bw1 = bwmorph(bw, 'fill');  figure; imshow(bw1);

SE = strel('disk',1,0);
bw1 = imfill(imdilate(bwmorph(bw, 'clean', inf), SE), 'holes');
figure; imshow(bw1);

[bw2 AreaTarget] = bwgetlargestblob(bw1);
figure; imshow(bw2);

%% extract posture and put it in animals' center of frame

figure; imshow(bw2); %this is already cropped

S = regionprops(bw2, {'Centroid', 'Orientation'});

figure; imshow(bw2); hold on; 
plot(S.Centroid(1), S.Centroid(2), 'r.');

%% analyze head and tail regions
