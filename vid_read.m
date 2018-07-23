% 7/22/2018 JJJ: created
function tmr = vid_read(vidobj, viF)
nThreads = 4; % number of parallel threads to run for loading video

fprintf('Loading video\n\t'); t1=tic;
tmr = zeros(vidobj.Height, vidobj.Width, numel(viF), 'uint8');

% parfor loading
if numel(viF)==1 || nThreads == 1
    fParfor = 0;
elseif median(diff(viF)) == 1
    fParfor = 0;
else
    fParfor = license('test', 'Distrib_Computing_Toolbox');
end
if fParfor
    fprintf('using parfor\n\t');
    try
        parfor (iF1=1:numel(viF), nThreads)
            img = read(vidobj, viF(iF1));
            tmr(:,:,iF1) = img(:,:,1);
            fprintf('.');
        end
    catch
        fprintf('parfor failed, retrying using for loop\n\t');        
        fParfor = 0;
    end %try
end
if ~fParfor
    for iF1=1:numel(viF)
        img = read(vidobj, viF(iF1));
        tmr(:,:,iF1) = img(:,:,1);
        fprintf('.');
    end
end
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func