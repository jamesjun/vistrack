function [vl_approach, vl_depart, vrRat] = getRegion(vs, mlMask, xy0)
vrX = poolVecFromStruct(vs, 'vrX');
vrY = poolVecFromStruct(vs, 'vrY');

vl = mlMask(sub2ind(size(mlMask), round(vrY), round(vrX)));

if nargin >= 3
    % approach is vl, depart is vl1
    vrR = sqrt((vrX - xy0(1)).^2 + (vrY - xy0(2)).^2);
    vrDR = differentiate3(vrR, .01);
    vrRat = vrDR ./ vrR;
    vlFast = abs(vrRat) > nanstd(vrRat(vl));
    vl_approach =  vl & (vrDR<0) & vlFast;
    vl_depart =  vl & (vrDR>0) & vlFast;
else
    vl_approach = vl;
    vl_depart = [];
end

if nargout == 0
%     figure; 
    hold on;
    
    for i=1:numel(vs)
        if iscell(vs)
            vs1 = vs{i};
        else
            vs1 = vs(i);
        end
        vrX = poolVecFromStruct(vs1, 'vrX');
        vrY = poolVecFromStruct(vs1, 'vrY');
        vl = mlMask(sub2ind(size(mlMask), round(vrY), round(vrX)));
        
        vrR = sqrt((vrX - xy0(1)).^2 + (vrY - xy0(2)).^2);
        vrDR = differentiate3(vrR);
        vl_approach =  vl & (vrDR<0);
        vl_depart =  vl & (vrDR>0);
    
        plot(vrX, vrY, 'k');
        plot(vrX(vl_approach), vrY(vl_approach), 'b.');
        plot(vrX(vl_depart), vrY(vl_depart), 'r.');
%         plot(vrX(end), vrY(end), 'r.');
%         plot(vrX(1), vrY(1), 'r.');
    end
end