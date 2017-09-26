function [ fname, dname, ext] = getFname( fullpath )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

istart = find(fullpath == '/', 1, 'last') + 1;
if isempty(istart)
    istart = find(fullpath == '\', 1, 'last') + 1;
end
if isempty(istart)
    istart = 1;
end

iend = find(fullpath == '.', 1, 'last') - 1;
if isempty(iend), iend = numel(fullpath); end

fname = fullpath(istart:iend);
if istart > 1
    dname = fullpath(1:istart-1);
else
    dname = '';
end

if iend < numel(fullpath)
    ext = fullpath(iend+1:end);
end

end

