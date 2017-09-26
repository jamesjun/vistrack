function disperr(vcMsg)
% ask user to email jrclust@vidriotech.com ? for the error ticket?
dbstack('-completenames'); 
if nargin==0
    vcMsg = lasterr();
end
fprintf(2, '%s\n', vcMsg);
try gpuDevice(1); disp('GPU device reset'); catch, end;
end