%GUI_ESACCADES

% SCRIPT_CUSTOM

% load the ESAC file and superimpose the location with the esaccades
% amplitudes.

LOADSETTINGS;
vrXh = handles.XC(:,2); %get the head location
vrYh = handles.YC(:,2); %get the head location
TC = handles.TC;

vrESAC = handles.ESAC.vr;
vxESAC = handles.ESAC.vx;
vyESAC = handles.ESAC.vy;
vtESAC = handles.ESAC.vt; 
viESAC = round(interp1(TC, 1:numel(TC), vtESAC, 'spline'));
threshESAC = quantile(vrESAC, .95);
viSel = find(vrESAC >= threshESAC);

figure; 
imshow(gray2rgb(handles.img0)); title('E-saccades locations');
hold on;

mrC = plotrainbow(vrXh, vrYh);

for i = 1:numel(viSel)
    iESAC = viSel(i);
    iTime = viESAC(iESAC);
    plot(vxESAC(iESAC), vyESAC(iESAC), 'o', 'color', mrC(iTime,:), 'MarkerSize', 8, 'LineWidth', 2);
end

return;

%location of ESAC

% prefix = getSpike2Prefix();
% sCh = getfield(ADC, sprintf('%s_Ch%d', prefix, ADC_CH));
% EODR = sCh.values;
% 
% for iframe=1:TRAJ_STEP:nframes
%     XI = interp1(2:nxy, XC(iframe,2:end), 2:.1:nxy, 'spline');
%     YI = interp1(2:nxy, YC(iframe,2:end), 2:.1:nxy, 'spline');
%     plot(XI, YI, 'color', mrColor(iframe,:));
%     plot(XI(1), YI(1), 'o', 'color', mrColor(iframe,:));
% end
