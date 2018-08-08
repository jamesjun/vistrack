% 7/22/2018 JJJ: created
function tmr = vid_read(vidobj, viF, nSkip_img)
if nargin<3, nSkip_img = []; end
if isempty(nSkip_img), nSkip_img = 1; end
if isempty(viF), viF = 1:vidobj.NumberOfFrames; end

nFrames_parfor = 300;
nThreads = 4; % number of parallel threads to run for loading video
    
fprintf('Loading video (%s: %d-%d, %d frames)\n', ...
    vidobj.Name, viF(1), viF(end), numel(viF)); 
t1=tic;
if nSkip_img == 1
    tmr = zeros(vidobj.Height, vidobj.Width, numel(viF), 'uint8');
else
    n1 = numel(1:nSkip_img:vidobj.Height);
    n2 = numel(1:nSkip_img:vidobj.Width);
    tmr = zeros(n1, n2, numel(viF), 'uint8');
end

% parfor loading
if numel(viF)<=nFrames_parfor || nThreads == 1
    fParfor = 0;
elseif median(diff(viF)) == 1
    fParfor = 0;
else
    fParfor = license('test', 'Distrib_Computing_Toolbox');
end
if fParfor
    fprintf('\tusing parfor\n');
    try
        parfor (iF1=1:numel(viF), nThreads)
            tmr(:,:,iF1) = read_(vidobj, viF(iF1), nSkip_img);
        end
    catch
        fprintf('parfor failed, retrying using for loop\n\t');        
        fParfor = 0;
    end %try
end
if ~fParfor
%     if all(diff(viF)==1)
%         tmr = read_(vidobj, viF([1,end]), nSkip_img);
% %         tmr = squeeze(tmr(:,:,1,:));        
%     else
        fprintf('\t');
        for iF1=1:numel(viF)
            tmr(:,:,iF1) = read_(vidobj, viF(iF1), nSkip_img);
%             tmr(:,:,iF1) = img(:,:,1);
%             fprintf('.');
        end
        fprintf('\n');
%     end
end
fprintf('\ttook %0.1fs\n', toc(t1));
end %func


function img = read_(vidobj, iFrame, nSkip_img)
if nargin<3, nSkip_img = 1; end

img = read(vidobj, iFrame);
img = img(:,:,1);
if nSkip_img>1
    img = binned_image_(img, nSkip_img, 0);
end
end %func


function img1 = binned_image_(img, nSkip, iMode)
% iMode: set to 0:averaging, 1:fast, 2:max
if nargin<3, iMode = 1; end
if ndims(img)==3, img = img(:,:,1); end
if iMode == 1 %fast mode
    img1 = img(1:nSkip:end, 1:nSkip:end); % faster
else
    dimm1 = floor(size(img)/nSkip);
    viY = (0:dimm1(1)-1) * nSkip;
    viX = (0:dimm1(2)-1) * nSkip;
    switch iMode 
        case 0
            img1 = zeros(dimm1, 'single');
            for ix = 1:nSkip
                for iy = 1:nSkip % TODO: reshape and use adjacent elements
                    img1 = img1 + single(img(viY+iy, viX+ix));
                end
            end
            img1 = img1 / (nSkip*nSkip);
            if isa(img, 'uint8'), img1 = uint8(img1); end
        case 2
            img1 = zeros(dimm1, 'like', img);
            for ix = 1:nSkip
                for iy = 1:nSkip
                    img1 = max(img1, img(viY+iy, viX+ix));
                end
            end
    end
end
end %func


