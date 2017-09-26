% plot ESAC
handles = load('F:\\2013.2S\\08-19_E02\\E02A1_Track');

vrXh = handles.XC(:,2);
vrYh = handles.YC(:,2);
vrESAC = ESAC.vr;
vxESAC = ESAC.vx;
vyESAC = ESAC.vy;

% find the location where E-saccades occur
threshESAC = quantile(vrESAC, .95);
viSel = find(vrESAC >= threshESAC);

figure; 
imshow(gray2rgb(handles.img0)); title('E-saccades locations');
hold on;

mrC = plotrainbow(vrXh, vrYh);

for i = 1:numel(viSel)
    iESAC = viSel(i);
    iTime = viESAC(iESAC);
    plot(vxESAC(iESAC), vyESAC(iESAC), 'o', 'color', mrC(iTime,:), 'MarkerSize', 8);
end

return;