classdef Trial
    %UNTITLED2 Summary of this class goes here
    %   XY already rotated and interpolated at 100 Frames/s

    properties (Constant)
        xyf = [789, 681];    rf = 1; %cm, radius
        xy1 = [966, 418];    r1 = 2.216*2.54/2; %*1.1222; %cm, radius
        xy2 = [975, 790];    r2 = 3.545*2.54/2; %*1.1222; %cm, radius
        xy3 = [604, 799];    r3 = 4*2.54/2; %cm, radius
        xy4 = [600, 428];    r4 = 3*2.54/2; %cm, radius
        xy0_ref = [787.0169, 605.6858]; %center
        pixpercm = 1053.28/(sqrt(2)*100);
    end

    
    properties
        TEOD, EODR, EODA, EODRq, EODRz, vtEOD, XH, YH, VEL, ACC, ANG, AVEL, ...
            HANG, HAVEL, viESAC, vrESAC, vtESAC, vxESAC, vyESAC, ...
            img0, dataID, pathLen_cm, xyFood, xy0, xyStart, vlZone, ...
            Rfood, Rcentre, vrX, vrY, mrX, mrY, iAnimal, duration, vcError;
    end
    
    
    methods
        function obj = Trial(S)
            if nargin == 0, return; end
            csFields = fieldnames(S);
            for iField = 1:numel(csFields)
                vcField = csFields{iField};
                try
                    eval(sprintf('field = S.%s;', vcField));
                    if size(field, 2) == 1 && size(field, 1) > 1, field=field'; end
                    eval(sprintf('obj.%s = field;', vcField, vcField));
                catch
                    disp(lasterr);
                end
            end
        end %Trial
        

        % Unit in cm
        function vrD = calcDistFrom(obj, vcFrom)
            vrX = obj.vrX - (obj.xy0(1) - Trial.xy0_ref(1));
            vrY = obj.vrY - (obj.xy0(2) - Trial.xy0_ref(2));
            
            switch lower(vcFrom)
                case 'food'
                    vrD = sqrt((vrX-Trial.xyf(1)).^2 + (vrY-Trial.xyf(2)).^2) / Trial.pixpercm - Trial.rf;
                case 'landmark1'
                    vrD = dist2square((vrX-Trial.xy1(1))/Trial.pixpercm, (vrY-Trial.xy1(2))/Trial.pixpercm, Trial.r1);
                case 'landmark2'
                    vrD = dist2square((vrX-Trial.xy2(1))/Trial.pixpercm, (vrY-Trial.xy2(2))/Trial.pixpercm, Trial.r2);
                case 'landmark3'
                    vrD = sqrt((vrX-Trial.xy3(1)).^2 + (vrY-Trial.xy3(2)).^2) / Trial.pixpercm - Tiral.r3;
                case 'landmark4'
                    vrD = sqrt((vrX-Trial.xy4(1)).^2 + (vrY-Trial.xy4(2)).^2) / Trial.pixpercm - Tiral.r4;
                case 'center'
                    vrD = sqrt((vrX-Trial.xy0_ref(1)).^2 + (vrY-Trial.xy0_ref(2)).^2) / Trial.pixpercm;
                otherwise
                    error('calcDistFrom-invalid location-%s', vcFrom);
            end
        end %calcDistFrom
    end %methods
end

