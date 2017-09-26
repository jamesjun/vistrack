function dimg = uint8_diff(img0c, img, fBright)
% dimg highlights dark by default and removes brighter feature. fBright rever5ses it.
if ndims(img0c)==3
    if ~fBright
%         dimg = uint8(ceil(mean(single(img0c) - single(img),3)));  
        dimg = (max(img0c - img, [], 3));
    else
%         dimg = uint8(ceil(mean(single(img) - single(img0c),3)));  
        dimg = (max(img - img0c, [], 3));
    end
else
    if ~fBright
        dimg = uint8(img0c - img);
    else
        dimg = uint8(img - img0c);
    end    
end
