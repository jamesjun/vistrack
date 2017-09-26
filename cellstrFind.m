function vi = cellstrFind(cellstr1, str1)
ci = strfind(cellstr1, str1);

vi = [];
for i=1:numel(ci)
    if ~isempty(ci{i})
        vi(end+1) = i;
    end
end


