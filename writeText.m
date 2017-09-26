function writeText(fname, csString)
fid = fopen(fname, 'w');
for i=1:numel(csString)-1
    fprintf(fid, '%s\n', csString{i});
end
fprintf(fid, '%s', csString{end});
fclose(fid);