function vrE = calcAngErrSign(vrA, vrB)
% rad

vrXa = cos(vrA(:));
vrYa = sin(vrA(:));
vrXb = cos(vrB(:));
vrYb = sin(vrB(:));
n = numel(vrA);

mrVa = [vrXa, vrYa, zeros(n,1)];
mrVb = [vrXb, vrYb, zeros(n,1)];
vrE = cross(mrVa, mrVb);
vrE = (asin(vrE(:,3)));