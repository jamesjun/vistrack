function data = differentiate3(data, dt)
%data = differentiate5(data, dt)
%data: timeseries to differentiate
%dt: time step (default of 1)
% http://en.wikipedia.org/wiki/Numerical_differentiation
dimm = size(data);
data=data(:)';
data = filter([1 0 -1], 2, data);
data = data(3:end);
data = [data(1), data, data(end)];

if nargin > 1
    data = data / dt;
end

if dimm(1) == 1 %row vector
    data=data(:)';
else
    data=data(:); %col vector
end