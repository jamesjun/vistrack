function img0=mask_copy(img1,img2,MASK)
img0 = img1;    
if ndims(img1)==2    
    if ~isempty(img2)
        img0(MASK) = img2(MASK);   
    else
        img0(MASK) = 0;
    end
else
    for i=1:size(img1,3)        
        img11 = img1(:,:,i);
        if ~isempty(img2)
            img21 = img2(:,:,i);
            img11(MASK) = img21(MASK);   
        else
            img11(MASK) = 0;
        end
        img0(:,:,i) = img11;
    end
end