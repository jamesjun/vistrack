function data = differentiate5(data, dt)
%data = differentiate5(data, dt)
%data: timeseries to differentiate
%dt: time step (default of 1)
% http://en.wikipedia.org/wiki/Numerical_differentiation
data=data(:)';
data = filter([-1 8  0 -8 1], 12, data);
data = data(5:end);
data = [data(1) * ones(1,2), data, data(end) * ones(1,2)];

if nargin > 1
    data = data / dt;
end