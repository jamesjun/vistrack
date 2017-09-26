function data = differentiate7(data, dt)
%data = differentiate5(data, dt)
%data: timeseries to differentiate
%dt: time step (default of 1)
% http://en.wikipedia.org/wiki/Numerical_differentiation
data=data(:)';
data = filter([1 -9 45 0 -45 9 -1], 60, data);
data = data(7:end);
data = [data(1) * ones(1,3), data, data(end) * ones(1,3)];

if nargin > 1
    data = data / dt;
end