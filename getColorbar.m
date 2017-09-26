function mrRGB = getColorbar(dimm)

mrRGB = zeros(dimm(1), dimm(2), 3);
vrColor = jet(dimm(1));

for i=1:dimm(1)
    j = dimm(1) - i + 1;
    mrRGB(i,:,1) =  vrColor(j,1);
    mrRGB(i,:,2) =  vrColor(j,2);
    mrRGB(i,:,3) =  vrColor(j,3);
end

