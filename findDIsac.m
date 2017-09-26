function vi = findDIsac(vrDI, thresh)
if nargin < 2
    thresh = [];
end

if isempty(thresh)
    thresh = -.04e-3;
end

vl = vrDI < thresh;
via = find(diff(vl) > 0) + 1;
vib = find(diff(vl) < 0) + 1;
vib(vib < via(1)) = []; % remove prev transitions
n = min(numel(via), numel(vib));
via = via(1:n);
vib = vib(1:n);

vnLen = vib - via;
vi = via(vnLen > 1);

if nargout == 0
    vrX = 1:numel(vrDI);
    figure; hold on;
    plot(vrX, vrDI, 'b.-');
    plot(vrX([1,end]), thresh*[1 1], 'k:');
    plot(vi, vrDI(vi), 'r.');
end
