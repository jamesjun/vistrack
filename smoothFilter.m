function X = smoothFilter( X, order)
%smoothFilter: based on smooth midpoints

    if nargin < 2
        order = 3;
    end

    order = round(order);
    A = pascal(order);
    coeff = zeros(order, 1);
    for i=1:order
        coeff(i) = A(i, order-i+1);
    end
    coeff = coeff / 2.^(order-1);


    X = filterSizeMatch(X, coeff);
end

