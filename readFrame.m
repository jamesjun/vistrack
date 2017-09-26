function mrImg = readFrame(vidobj, iframe)
CACHESIZE = 200;
persistent FLIM trImg;

fRebuild = 0;
if isempty(FLIM)
    fRebuild = 1;
else    
    if iframe > FLIM(2) || iframe < FLIM(1)
        fRebuild = 1;
    end
end

if fRebuild    
    FLIM = iframe + [0 CACHESIZE];
    FLIM(2) = min(FLIM(2), vidobj.NumberOfFrames);
    fprintf('loading video frames %d-%d\n', FLIM(1), FLIM(2));
    trImg = read(vidobj, FLIM);
end

mrImg = trImg(:,:,:, iframe - FLIM(1) + 1);
mrImg = mrImg(:,:,1);