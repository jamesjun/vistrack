function ESAC = calcESAC(handles)

% load the ESAC file and superimpose the location with the esaccades
% amplitudes.

LOADSETTINGS;
vrXh = handles.XC(:,2); %get the head location
vrYh = handles.YC(:,2); %get the head location
TC = handles.TC;

%% load ADC file (EODA)
[dirname, filename, ~] = fileparts(handles.vidFname);
fname_esac = sprintf('%s\\%s_EODA.mat', dirname, filename);
try
    ADC_ESAC = load(fname_esac);
catch
    disp('EODA not saved.');
    ESAC = [];
    return;
end
prefix = getSpike2Prefix(ADC_ESAC);
sCh = getfield(ADC_ESAC, sprintf('%s_Ch%d', prefix, ADC_CH_ESAC));
vrESAC = sCh.values;
vtESAC = sCh.times;
viSel = vtESAC >= TC(1) & vtESAC < TC(end);
vtESAC = vtESAC(viSel);
vrESAC = vrESAC(viSel);

% figure; plot(vtESAC, vrESAC, '.');
% figure; ksdensity(log10(vrESAC));

vxESAC = interp1(TC, vrXh, vtESAC, 'spline');
vyESAC = interp1(TC, vrYh, vtESAC, 'spline');

ESAC.vt = vtESAC;
ESAC.vr = vrESAC;
ESAC.vx = vxESAC;
ESAC.vy = vyESAC;

%% plot
if nargout > 0
    return;
end

% find the location where E-saccades occur

threshESAC = quantile(vrESAC, .95);
viSel = find(vrESAC >= threshESAC);

figure; 
imshow(gray2rgb(handles.img0)); title('E-saccades locations');
hold on;

plotrainbow(vrXh, vrYh);

for i = 1:numel(viSel)
    iESAC = viSel(i);
    plot(vxESAC(iESAC), vyESAC(iESAC), 'wo');
end

return;

