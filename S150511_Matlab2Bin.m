% tritrode to klustakwik .bin format
vcDir = '/Users/junj10/Dropbox (HHMI)/MladenTritrodePatchData/01092013/';
[mrWav, S] = importAbf(vcDir, 'fFibo', 0, 'fSpont', 1, 'timeLim', [3 inf], 'fPlot', 0);
S.vcFname

%%
if iscell(S)
    S1 = S{1}; 
else
    S1 = S;
end

nSamples = size(mrWav,1);
fid = fopen('tritrode3.dat', 'w');
step = 0.006103515761424 / 10;
fwrite(fid, int16(toVec(mrWav(:, S1.chExt)'/step)), 'int16', 'l');
fclose(fid);
