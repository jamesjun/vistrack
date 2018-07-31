% 7/22/2018 JJJ: created
function tmr = vid_read(vidobj, viF)
nThreads = 4; % number of parallel threads to run for loading video

fprintf('Loading video (%s: %d-%d, %d frames)\n', ...
    vidobj.Name, viF(1), viF(end), numel(viF)); 
t1=tic;
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
    fprintf('\tusing parfor\n');
    try
        parfor (iF1=1:numel(viF), nThreads)
            img = read(vidobj, viF(iF1));
            tmr(:,:,iF1) = img(:,:,1);
%             fprintf('.');
        end
    catch
        fprintf('parfor failed, retrying using for loop\n\t');        
        fParfor = 0;
    end %try
end
if ~fParfor
    if all(diff(viF)==1)
        tmr = read(vidobj, viF([1,end]));
        tmr = squeeze(tmr(:,:,1,:));        
    else
        fprintf('\t');
        for iF1=1:numel(viF)
            img = read(vidobj, viF(iF1));
            tmr(:,:,iF1) = img(:,:,1);
            fprintf('.');
        end
        fprintf('\n');
    end
end
fprintf('\ttook %0.1fs\n', toc(t1));
end %func