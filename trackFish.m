function [XC, YC, AC, Area, S, MOV, XC_off, YC_off] = trackFish(S, FLIM)
% S.{AreaTarget, ThreshLim, img0, fShow, vec0, thresh, WINPOS}
% AC: degree unit
% XC, YC: pixel unit

MEMLIM = 300; %number of frames to load to memory at a time
nRetry = 3; %number of retries for loading a video file

% Parse input variables
WINPOS = S.WINPOS;
vecPrev = S.vec0;
thresh = S.thresh;
nframes = diff(FLIM) + 1;

% Allocate output arrays
XC = nan(nframes,6);
YC = nan(nframes,6);
AC = nan(nframes,5);
Area = nan(nframes,1);
height = diff(WINPOS([3 4])) + 1;
width = diff(WINPOS([1 2])) + 1;
MOV = zeros(height, width, nframes, 'uint8');
XC_off = nan(nframes, 1);
YC_off = nan(nframes, 1);

if nargin < 2
    FLIM = [1 S.vidobj.NumberOfFrames];
end

% Call itself recursively if the number of frames is over the memory limit
if nframes > MEMLIM                
    for iF=FLIM(1):MEMLIM:FLIM(2)
        FLIM1 = [iF, iF+MEMLIM-1];
        FLIM1(2) = min(FLIM1(2), FLIM(2));
        LLIM1 = FLIM1 - FLIM(1) + 1;
        L = LLIM1(1):LLIM1(2);
        try
            [XC(L,:), YC(L,:), AC(L,:), Area(L,:), S1, MOV(:,:,L), XC_off(L), YC_off(L)] ...
                = trackFish(S, FLIM1);
        catch
            disperr();
        end
        S = S1;
    end
    return;
end

tic; %start the timer

% Load video frames to the memory
for itry=1:nRetry
    try
        h=msgbox(sprintf('Tracking frames %d ~ %d... (close to cancel)', FLIM(1), FLIM(2)));
        IMG = read(S.vidobj, FLIM);
%         IMG = IMG(:,:,1,:); %use red channel only
        try close(h); 
        catch, error('Cancelled by user'); end;        
        break;
    catch
        disp(lasterr);
        fprintf('failed to load %d times. reloading...\n', itry);
        S.vidobj = VideoReader(S.vidFname);
    end
end
if itry == nRetry
    error('video load failure');
end

if S.fShow
    hfig = figure;
end

% Process each frame
[hAx, hImg, hPlot1, hPlot2] = deal([]);
BWprev=[]; BWprev1=[];
for iF=1:nframes
    [img, dimg, dimg1] = getFrame_(IMG, iF, WINPOS, S.img0);
    BW0 = (dimg(:,:,1) > thresh);
    BW = imdilate(bwmorph(BW0, 'clean', inf), S.SE);
    BW1 = BW;
    iF1 = FLIM(1) + iF - S.FLIM(1);
    regions = regionprops(BW, {'Area', 'Centroid', 'Orientation', 'FilledImage', 'BoundingBox'});
    
    %remove blobs not touching the fish's blob from the prev. frame
    if numel(regions) > 1
        L = bwlabel(BW, 8);
        if iF1<=10 % first frame, pick closest to the center
            [iRegion] = region_nearest(regions, S.xy0 - WINPOS([1,3]));
            regions = regions(iRegion);
            BW = L==iRegion;
        else
            [iRegion] = region_largest(regions);
            regions = regions(iRegion);
            BW = L==iRegion;
        end
    end
    
    %Isolate the largest blob
    if isempty(regions)    
        BW=BWprev; % revert to the 
        regions = regionprops(BW, {'Area', 'Centroid', 'Orientation', 'FilledImage', 'BoundingBox'});
        stats = largestBlob_(regions);
    else
        stats = regions;
    end
    
    ang = -stats.Orientation;
    area = stats.Area;
    xy_cm = stats.Centroid;    
    %Check for the orientation flip
    vec = [cos(deg2rad(ang)), sin(deg2rad(ang))];
    if dot(vec, vecPrev) < 0
        ang = mod(ang + 180, 360);
        if ang>180, ang=ang-360; end
        vec = [cos(deg2rad(ang)), sin(deg2rad(ang))];
        stats.Orientation = -ang;
    end    
    
    %Compute posture
    try
        [XY, ANG, xy_names, ang_names] = blobPosture(stats);
    catch
        img1 = img;
        img1(bwperim(BW)) = 255;
        figure(101); imshow(img1); title('Blob processing error'); drawnow;
        fprintf(2, 'Blob processing error\n'); continue;
    end
    
    %Display output
    if S.fShow
        img1 = img; 
        img1(bwperim(BW)) = 255;
          
        try
            delete(hPlot1);
            delete(hPlot2);
        catch
        end
        if isempty(hImg)
            figure(hfig);
            hImg = imshow(img1);
            hAx = gca;
        else
            set(hImg,'CData', img1);
        end
        hold on;
        hPlot1 = plot(hAx, XY(2, 1), XY(2, 2), 'go', ...
            XY(3:end, 1), XY(3:end, 2), 'mo', 'LineWidth', 2);
        
        % interpolated curve
        nxy = size(XY,1);
        X1 = interp1(2:nxy, XY(2:end, 1), 2:.1:nxy, 'spline');
        Y1 = interp1(2:nxy, XY(2:end, 2), 2:.1:nxy, 'spline');
        hPlot2 = plot(hAx, X1, Y1, 'm-', XY(1,1), XY(1,2), 'g+', 'LineWidth', 2); %Mark the centroid
        title(hAx, sprintf('Frame = %d', iF + FLIM(1) - 1));
        drawnow;
    end
        
    %Save output
    Area(iF) = area;
    xy_off = [WINPOS(1), WINPOS(3)] - [1, 1];
    XC(iF,:) = round(XY(:,1)' + xy_off(1));
    YC(iF,:) = round(XY(:,2)' + xy_off(2));
    AC(iF,:) = normAng(ANG');
    MOV(:,:,iF) = img;
    XC_off(iF) = xy_off(1); 
    YC_off(iF) = xy_off(2);
    
    %Update the bounding box
%     xy_center = xy_cm + xy_off;
    xy_center = getMedian(XC(:,1), YC(:,1), iF); %fault tolerent
    [WINPOS, ~] = getBoundingBoxPos(xy_center, size(S.img0), size(BW));    
    
    %adjust intensity threshold
    if area > S.AreaTarget*1.1
        thresh = min(thresh+1, S.ThreshLim(2));
    elseif area < S.AreaTarget*.9
        thresh = max(thresh-1, S.ThreshLim(1));
    end
    
    %next orientation vector
%     vecPrev = vec;
    vecPrev = getMedian(cos(deg2rad(AC(:,1))), sin(deg2rad(AC(:,1))), iF);
    vecPrev = vecPrev ./ norm(vecPrev);

    BWprev1 = BWprev;
    BWprev = BW;
    if isempty(BWprev1), BWprev1=BWprev; end
end %for

%Return the last iteration info
S.thresh = thresh;
S.vec0 = vec;
S.WINPOS = WINPOS;
S.xy_names = xy_names;
S.ang_names = ang_names;

%Measure the processing time
tdur = toc;
fprintf('Processed %0.1f images/sec, %s, Frames: [%d ~ %d]\n', ...
        nframes/tdur, S.vidobj.Name, FLIM(1), FLIM(2));
    
try close(hfig); catch, end;

end %func


%--------------------------------------------------------------------------
function xy = getMedian(vrX, vrY, idx)
n = 5;

idxrng = [idx-n+1:idx];
idxrng(1) = max(idxrng(1), 1);
try
    xy = [median(vrX(idxrng)), median(vrY(idxrng))];
catch
    xy = [vrX(idx), vrY(idx)];
end
end

%--------------------------------------------------------------------------
function [img, dimg, dimg1] = getFrame_(IMG, iF, WINPOS, img0)
%GETFRAME Get image frame and crop and subtract the background
% [img] = getFrame(IMG, iF) %Obtain current frame from array of images
% [img] = getFrame(IMG, iF, WINPOS) %crops the image
% [img, dimg] = getFrame(IMG, iF, WINPOS, img0) %crop and subtract
%  background

if nargin < 3
    WINPOS = [1, size(IMG, 2), 1, size(IMG,1)];
end

img = IMG(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2), :, iF);

% Calculate the intensity difference (background subtraction)
if nargout >= 2
    img0c = img0(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2), :);
    dimg = uint8_diff(img0c, img, 0);
    if nargout>=3
        dimg1 = uint8_diff(img0c, img, 1);
    end
%    dimg = img0(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2)) - img;
end
img = uint8(single(mean(img,3)));
end %func


%--------------------------------------------------------------------------
%NORMANG Normalize the angle to range between -180 to 180 degrees
% [A] = normAng(A)
function A = normAng(A)
A = mod(A, 360);
A(A>180) = A(A>180) - 360;
end %func


%--------------------------------------------------------------------------
%LARGESTBLOB Return stats for the largest blob
% [stat, idx, area] = largestBlob(stats)
function [stat, idx, area] = largestBlob_(stats)
[area, idx] = max([stats.Area]);
stat = stats(idx);
end

%--------------------------------------------------------------------------
% Measure posture angle from a binary blob 
function [XY, ANG, xy_names, ang_names] = blobPosture(stats)

BW0 = stats.FilledImage;
ang0 = -stats.Orientation; %counter-clockwise is positive in pixel 
xy_ref = [stats.BoundingBox(1), stats.BoundingBox(2)];
xy_ref = round(xy_ref + [stats.BoundingBox(3), stats.BoundingBox(4)]/2);
xy_cm = stats.Centroid;

%Rotate original image (BW0) parallel to the major axis
% BWr = imrotate(BW0, ang0);
BWr = imdilate(imrotate(BW0, ang0), ones(3));
xy_r = round([size(BWr, 2), size(BWr, 1)]/2); %rotation center. use this as a reference
% stats_r = regionprops(BWr, 'Area', 'BoundingBox'); %have area for safety
% [~, idx] = max([stats_r.Area]); 
% stats_r = stats_r(idx);
stats_r = largestBlob_(regionprops(BWr, 'BoundingBox', 'Area'));

% find head CM
BWh = BWr(1:end, xy_r(1):end);
% stats_h = largestBlob(regionprops(BWh, 'Orientation', 'Centroid', 'Image', 'Area'));
stats_h = largestBlob_(regionprops(BWh, 'Orientation', 'Centroid', 'Area'));
ang_h = -stats_h.Orientation;
xy_hm = round(stats_h.Centroid);
xy_hm(1) = xy_hm(1) + xy_r(1) - 1;
xy_hm(2) = median(find(BWr(:, xy_hm(1))));

% find tail CM
BWt = BWr(1:end, 1:xy_r(1));
% stats_t = largestBlob(regionprops(BWt, 'Orientation', 'Centroid', 'Image', 'Area'));
stats_t = largestBlob_(regionprops(BWt, 'Orientation', 'Centroid', 'Area'));
ang_t = -stats_t.Orientation;
xy_tm = round(stats_t.Centroid);
xy_tm(2) = median(find(BWr(:, xy_tm(1))));

%find middle points
xy_m = [xy_r(1), nan];
xy_m(2) = median(find(BWr(:, xy_r(1))));

%find tail tip
xy_t = [ceil(stats_r.BoundingBox(1)) , nan];
% xy_t(2) = find(BWr(:,xy_t(1)), 1, 'first');
xy_t(2) = median(find(BWr(:, xy_t(1))));

%find head tip
% xy_h = [floor(sum(stats_r.BoundingBox([1, 3]))) , nan];
% dx = xy_h(1) - xy_hm(1);
% xy_h(2) = round(dx * tan(deg2rad(ang_h)) + xy_hm(2));
% xy_h(2) = round(find(BWr(:,xy_h(1)), 1, 'first'));

%find head tip
BWrr = imrotate(BWr, ang_h);
stats_rr = largestBlob_(regionprops(BWrr, 'Orientation', 'BoundingBox', 'Area'));
xy_rr = [size(BWrr, 2), size(BWrr, 1)]/2; %rotation center. use this as a reference

xy_h = [floor(sum(stats_rr.BoundingBox([1, 3]))) , nan];
try
%     xy_h(2) = median(find(BWrr(:, xy_hm(1))));
    xy_h(2) = median(find(BWrr(:, xy_h(1))));
catch
    disperr();
end
% xy_h(2) = round(median(find(BWrr(:, xy_h(1)))));

% compute angles
vec1 = xy_hm - xy_m;
vec2 = xy_m - xy_t;
ang_tb = rad2deg(atan2(vec2(2),vec2(1)) - atan2(vec1(2),vec1(1)));

%format output
ang_names = {'CoM', 'head-mid', 'tail-mid', 'body-bend', 'tail-bend'};
ANG = zeros(numel(ang_names), 1);
ANG(1) = ang0;
ANG(2) = ang_h + ang0;
ANG(3) = ang_t + ang0;
ANG(4) = ang_t - ang_h;
ANG(5) = ang_tb;

%compute positions
xy_r = [size(BWr, 2), size(BWr, 1)]/2; %do not round for higher precision
xy_names = {'CoM', 'head', 'head-mid', 'mid', 'tail-mid', 'tail'};
XY = zeros(numel(xy_names), 2);
XY(1,:) = xy_cm;
XY(2,:) = xy_ref + rotatexy(xy_h - xy_rr, ang0 + ang_h)';
XY(3,:) = xy_ref + rotatexy(xy_hm - xy_r, ang0)';
XY(4,:) = xy_ref + rotatexy(xy_m - xy_r, ang0)';
XY(5,:) = xy_ref + rotatexy(xy_tm - xy_r, ang0)';
XY(6,:) = xy_ref + rotatexy(xy_t - xy_r, ang0)';

end %func


%--------------------------------------------------------------------------
function [ xyp ] = rotatexy( xy, ang )
%ROTATEXY rotate a vector with respect to the origin, ang in degree
xy = xy(:);
CosA = cos(deg2rad(ang));
SinA = sin(deg2rad(ang));
M = [CosA, -SinA; SinA, CosA];
xyp = M * xy;
end


%--------------------------------------------------------------------------
function [ rad ] = deg2rad( deg )
%DEG2RAD convert an angle from degrees to radians
rad = deg / 180 * pi;
end

%--------------------------------------------------------------------------
function [ deg ] = rad2deg( rad )
%RAD2DEG convert an angle from radians to degrees
deg = rad / pi * 180;
end
