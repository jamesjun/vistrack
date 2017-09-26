% visual tracking. adapted from Ms_inverseProblem
fShow = 1; %show display

% new 8chan HDvid
vidFname = 'D:\Expr\121212_spontvid8\121212_spontvid8ch000-1.wmv';
outFname = 'VISTRACK_121212_000.mat';  
CAMRES = 'HD';
TLIM0 = [10 5580];

switch CAMRES
    case 'SD'   %640x480 @20FPS
    winlen = 128*1.5;
    FPS_cam = 20.066995768688294;
    SE = strel('disk',1,0); %use disk size 1 for 640 480
    ThreshLim = [0 20];
    IM_THRESH = 4;
    IM_CONTRAST = .85;

    case 'HD'   %1600x1200 @15FPS
    winlen = 128*2.5;
    FPS_cam = 15.001551480959098;
    SE = strel('disk',3,0); %use disk size 1 for 640 480
    ThreshLim = [0 15];
    IM_THRESH = 5;
    IM_CONTRAST = .85;
end

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

%%
% intrng = stretchlim(img1, [0 .85]);
% img1s = imadjust(img1, intrng);
% img2s = imadjust(img2, intrng);
% [img00, MASK, xy_init, vec0, xy0] = make_BG_JJJ(img1s, img2s);

%% Preview initial condition
IM_CONTRAST = .9;
IM_THRESH = 5;

img0 = img00* IM_CONTRAST; %make it darker
thresh = IM_THRESH;

% Initial
[WINPOS, ~] = getBoundingBoxPos(xy_init, size(img0), winlen*[1 1]);
img = img1(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2));
dimg = uint8(img0(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2))-img); 
BW = imdilate(bwmorph((dimg > thresh), 'clean', inf), SE);
BW = imfill(BW, 'holes');
[BW, AreaTarget] = bwgetlargestblob(BW);


figure; 
subplot 231; hfig1 = imshow(img);  
hold on;    hax1 = gca; 
title('original');

subplot 232; hfig2 = imshow(dimg); 
hold on;    hax2 = gca; 
title('difference');

subplot 233; hfig3 = imshow(BW);
hold on;    hax3 = gca; 
title(sprintf('Area=%d, threshold=%d', AreaTarget, thresh));    

subplot 234;  
img4 = img; img4(~BW)=img4(~BW)*.5;     
hfig4 = imshow(img4); 
hold on; hax4 = gca; 
title('superimposed');
colormap gray;

subplot 235; hfig5 = imshow(img0);  
hold on;    hax5 = gca; 
title('Background image');
drawnow;

ans = questdlg('OK to proceed?', 'proceed?', 'yes', 'no', 'yes');
if strcmp(ans, 'no'), return;   end


%% track
%no smoothing. each frame
S = makeStruct( AreaTarget, ThreshLim, vec0, img0, thresh, WINPOS, SE, aviobj, vidFname, TLIM0, FLIM0, FLIM, fShow, xy0);
[XC, YC, AC, Area, S1] = trackFish(S, FLIM);

% save result of visual tracking
xy_names = S1.xy_names;
ang_names = S1.ang_names;
eval(sprintf('save %s TC XC YC AC Area S xy_names ang_names;', outFname));

if ~fShow
    return;
end

%% plot the trajectories

% filter and differentiate
%use the head-mid point as an x,y,a marker, and tail bending angle for body
%bending

nf = 3;
XCf_h = filtPos(XC(:,2),nf);    YCf_h = filtPos(YC(:,2),nf);
XCf_hm = filtPos(XC(:,3),nf);   YCf_hm = filtPos(YC(:,3),nf);
XCf_m = filtPos(XC(:,4),nf);    YCf_m = filtPos(YC(:,4),nf);
XCf_tm = filtPos(XC(:,5),nf);    YCf_tm = filtPos(YC(:,5),nf);
XCf_t = filtPos(XC(:,6),nf);    YCf_t = filtPos(YC(:,6),nf);

% figure; plot(TC,XC(:,3),'.:', TC,YC(:,3),'.:', TC,AC(:,2),'.:', TC,AC(:,5),'.:');
% hold on;
% plot(TC,XCf, TC,YCf, TC,ACf, TC,BCf);

% verify the tracking with movie
figure; fig = imshow(img0);
intlim = stretchlim(img00(MASK), [0 .9]);
for i=1:15:numel(TC)
    img = read(aviobj, FLIM(1) + i - 1);
    img = imadjust(img, intlim);
    set(fig, 'cdata', img(:,:,1)); hold on;
%     h1=plot(XCf(i) + cos(deg2rad(ACf(i))) * [50 0], ...
%          YCf(i) + sin(deg2rad(ACf(i))) * [50 0], 'b.-');
%     h2=plot(XCf(i) - cos(deg2rad(ACf(i)+BCf(i))) * [100 0], ...
%          YCf(i) - sin(deg2rad(ACf(i)+BCf(i))) * [100 0], 'r.-');     
    h1 = plot([XCf_h(i), XCf_m(i)], [YCf_h(i), YCf_m(i)], 'r.-');
    h2 = plot([XCf_m(i), XCf_t(i)], [YCf_m(i), YCf_t(i)], 'b.-');
    h3 = plot(XCf_m(i), YCf_m(i), 'g*');
    h4 = plot(XCf_hm(i), YCf_hm(i), 'r*');
    h5 = plot(XCf_tm(i), YCf_tm(i), 'b*');
    title(sprintf('t=%0.1f sec', TC(i)));
     drawnow;
     
     delete([h1 h2 h3 h4 h5]);
end
