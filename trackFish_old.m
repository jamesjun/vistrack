function [XC, YC, AC, Area, S] = trackFish(S, FLIM)
% S.{AreaTarget, ThreshLim, img0, fShow, vec0, thresh, WINPOS}
MEMLIM = 200; %number of frames to load to memory at a time
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
        [XC(L,:), YC(L,:), AC(L,:), Area(L,:), S] = trackFish(S, FLIM1);
    end
    return;
end

tic; %start the timer

% Load video frames to the memory
for itry=1:nRetry
    try
        h=msgbox(sprintf('Tracking frames %d ~ %d... (this will close automatically)', FLIM(1), FLIM(2)));
        IMG = read(S.vidobj, FLIM);
        try close(h); catch, end;
        IMG = IMG(:,:,1,:); %use red channel only
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
for iF=1:nframes
    [img, dimg] = getFrame(IMG, iF, WINPOS, S.img0);
    BW0 = (dimg > thresh);
    BW = imdilate(bwmorph(BW0, 'clean', inf), S.SE);
    BW = imclearborder(BW,8); %remove boundary touching    

    %Isolate the largest blob
    stats = largestBlob(...
            regionprops(BW, {'Area', 'Centroid', 'Orientation', 'FilledImage', 'BoundingBox'}));
    try
        ang = -stats.Orientation;    
        area = stats.Area;
        xy_cm = stats.Centroid;
    catch
       disp(lasterr);
    end
    
    %Check for the orientation flip
    vec = [cos(deg2rad(ang)), sin(deg2rad(ang))];
    if dot(vec, vecPrev) < 0
        ang = mod(ang + 180, 360);
        if ang>180, ang=ang-360; end
        vec = [cos(deg2rad(ang)), sin(deg2rad(ang))];
        stats.Orientation = -ang;
    end    
    
    %Compute posture
    [XY, ANG, xy_names, ang_names] = blobPosture(stats);
    
    %Display output
    if S.fShow
        img1 = img; 
        BW1 = bwperim(bwgetlargestblob(BW));
        img1(BW1) = 255;
                
        clf(hfig);
        figure(hfig);
        
        imshow(img1); hold on;
        plot(XY(2, 1), XY(2, 2), 'wo'); %head
        plot(XY(3:end, 1), XY(3:end, 2), 'ro');
        
        % interpolated curve
        nxy = size(XY,1);
        X1 = interp1(2:nxy, XY(2:end, 1), 2:.1:nxy, 'spline');
        Y1 = interp1(2:nxy, XY(2:end, 2), 2:.1:nxy, 'spline');
        plot(X1, Y1, 'r-');
        plot(XY(1,1), XY(1,2), 'g+'); %Mark the centroid
        title(sprintf('Frame = %d', iF + FLIM(1) - 1));
        
        drawnow;
    end
        
    %Save output
    Area(iF) = area;
    xy_off = [WINPOS(1), WINPOS(3)] - [1, 1];
    XC(iF,:) = round(XY(:,1)' + xy_off(1));
    YC(iF,:) = round(XY(:,2)' + xy_off(2));
    AC(iF,:) = normAng(ANG');
    
    %Update the bounding box
    [WINPOS, ~] = getBoundingBoxPos(xy_cm + xy_off, size(S.img0), size(BW));    
    
    %adjust intensity threshold
    if area > S.AreaTarget*1.1
        thresh = min(thresh+1, S.ThreshLim(2));
    elseif area < S.AreaTarget*.9
        thresh = max(thresh-1, S.ThreshLim(1));
    end  
    vecPrev = vec; %next orientation vector    
end %for

%Return the last iteration info
S.thresh = thresh;
S.vec0 = vec;
S.WINPOS = WINPOS;
S.xy_names = xy_names;
S.ang_names = ang_names;

%Measure the processing time
tdur = toc;
fprintf('Took %0.1f images/sec, %s, Frames: [%d ~ %d]\n', ...
        nframes/tdur, S.vidobj.Name, FLIM(1), FLIM(2));
    
try close(hfig); catch, end;

end %func


%--------------------------------------------------------------------------
function [img, dimg] = getFrame(IMG, iF, WINPOS, img0)
%GETFRAME Get image frame and crop and subtract the background
% [img] = getFrame(IMG, iF) %Obtain current frame from array of images
% [img] = getFrame(IMG, iF, WINPOS) %crops the image
% [img, dimg] = getFrame(IMG, iF, WINPOS, img0) %crop and subtract
%  background

if nargin < 3
    WINPOS = [1, size(IMG, 2), 1, size(IMG,1)];
end

img = IMG(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2), 1, iF);

% Calculate the intensity difference (background subtraction)
if nargin >= 4
   dimg = img0(WINPOS(3):WINPOS(4), WINPOS(1):WINPOS(2)) - img;
end

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
function [stat, idx, area] = largestBlob(stats)
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
BWr = imrotate(BW0, ang0);
xy_r = round([size(BWr, 2), size(BWr, 1)]/2); %rotation center. use this as a reference
stats_r = regionprops(BWr, 'Area', 'BoundingBox'); %have area for safety
[~, idx] = max([stats_r.Area]); 
stats_r = stats_r(idx);

% find head CM
BWh = BWr(1:end, xy_r(1):end);
stats_h = largestBlob(regionprops(BWh, 'Orientation', 'Centroid', 'Image', 'Area'));
ang_h = -stats_h.Orientation;
xy_hm = round(stats_h.Centroid);
xy_hm(1) = xy_hm(1) + xy_r(1) - 1;
xy_hm(2) = median(find(BWr(:, xy_hm(1))));

% find tail CM
BWt = BWr(1:end, 1:xy_r(1));
stats_t = largestBlob(regionprops(BWt, 'Orientation', 'Centroid', 'Image', 'Area'));
ang_t = -stats_t.Orientation;
xy_tm = round(stats_t.Centroid);
xy_tm(2) = median(find(BWr(:, xy_tm(1))));

%find middle points
xy_m = [xy_r(1), nan];
xy_m(2) = median(find(BWr(:, xy_r(1))));

%find tail tip
xy_t = [ceil(stats_r.BoundingBox(1)) , nan];
xy_t(2) = find(BWr(:,xy_t(1)), 1, 'first');

%find head tip
% xy_h = [floor(sum(stats_r.BoundingBox([1, 3]))) , nan];
% dx = xy_h(1) - xy_hm(1);
% xy_h(2) = round(dx * tan(deg2rad(ang_h)) + xy_hm(2));
% xy_h(2) = round(find(BWr(:,xy_h(1)), 1, 'first'));

%find head tip
BWrr = imrotate(BWr, ang_h);
stats_rr = largestBlob(regionprops(BWrr, 'Orientation', 'BoundingBox', 'Area'));
xy_rr = [size(BWrr, 2), size(BWrr, 1)]/2; %rotation center. use this as a reference

xy_h = [floor(sum(stats_rr.BoundingBox([1, 3]))) , nan];
xy_h(2) = median(find(BWrr(:, xy_hm(1))));
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