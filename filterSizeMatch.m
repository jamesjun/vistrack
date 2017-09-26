function [ X n1 n2] = filterSizeMatch( X, KERN)
%FILTERSIZEMATCH returns a filtered data which matches the size
%n1, n2: from and to where filter is applied

    nshift = round(numel(KERN)/2);
    Xfirst = X(1:nshift); 
    Xlast = X(end-nshift+1:end);  
    X = conv(X, KERN, 'same');
    X(1:nshift) = Xfirst;
    X(end-nshift+1:end) = Xlast;
    
    n1 = nshift+1;
    n2 = numel(X) - nshift;
end




