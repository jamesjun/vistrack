function [ y ] = filtPos( x, n, dimm)
%FILTPOS apply smoothing filter to the position (x,y)
%James J. Jun. Malerlab. 2012 Nov
% Median filter then averaging
%
%[Variables]
%x: input position
%n: filter window
%y: smoothed position

if nargin < 2, n = 15; end
if nargin < 3
    dimm = [];
end

%filter each dimm
if min(size(x)) > 1
    if dimm == 2
        x = x';
    end
    y = zeros(size(x));
    for iCol=1:size(x,2)
        y(:,iCol) = filtPos(x(:,iCol), n);
    end
    if dimm == 2
        y = y';
    end    
    return;
end

%make the median filter length odd
n1 = round(n/2); if mod(n1,2)==0, n1 = n1 + 1; end 

%deal with nan
IDX = find(isnan(x(2:end-1)))+1;
if ~isempty(IDX)
    for i=1:numel(IDX)
        idx = IDX(i);
        x(idx) = nanmean(x(idx-1:idx+1));
    end
end

y1 = medfilt1(x,n1);
y = filtfilt(ones(1,n), n, y1);
y(1:n) = nanmedian(x(1:n)); %fix the start
end

