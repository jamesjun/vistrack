function [y, y_det] = sortnorm(x)
[~, ix] = sort(x);
if ismatrix(ix)
    y = zeros(size(x));
    y1 = (1:size(x,1))' / size(x,1);
    for iC=1:size(x,2)
        y(ix(:,iC),iC) = y1;
    end
%     y(ix) = repmat(y1, [1, size(x,2)]);
    if nargout>=2
        y_det = bsxfun(@minus, y, y1);
    end    
else
    y1 = (1:numel(x))/numel(x);
    y(ix) = y1;
    if nargout>=2
        y_det = y - y1;
    end
end
if isa(x,'uint8'), y = uint8(y*255); end