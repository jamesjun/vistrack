 
% fname = 'D:\Expr\121212_spontvid8\121212_spontvid8ch000-1.wmv';   
% FLIM = [1 15*21];

fname = '121126_SPONT000-1.wmv';
FLIM = [4300 4510];

vobj = VideoReader(fname);
vmat = read(vobj, FLIM);
writerObj = VideoWriter('sample.avi', 'Uncompressed AVI');
writerObj.FrameRate = 15;
% writerObj.Quality = 100;
open(writerObj);
writeVideo(writerObj,vmat);
close(writerObj);
