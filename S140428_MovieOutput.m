fname = 'C:\expr\Raw data\E13A4p.wmv';
vid = VideoReader(fname);
MOV_FILEOUT = 'output.mp4';


% load data
LOADSETTINGS;
TC = handles.TC;
XC = handles.XC;
YC = handles.YC;
try
    ADC = handles.ADC;
    [EODR, TEOD, chName] = getSpike2Chan(ADC, ADC_CH_PLOT);
    fprintf('Loaded Spike2 %s (Ch%d)\n', chName, ADC_CH_PLOT);
    
    % Resample EOD Rate to camera time    
    RC = interp1(TEOD, EODR, TC);
%     RC = filtfilt(ones(1,TRAJ_STEP), TRAJ_STEP, RC);
catch
    disp(lasterr)
    errordlg('Load ADC file');
    return;
end
% figure; plot(TEOD, EODR, 'r.', TC, RC, 'b-');

% Plot the EOD color representation
writerObj = VideoWriter(MOV_FILEOUT, 'MPEG-4');
set(writerObj, 'FrameRate', handles.FPS); %30x realtime
% set(writerObj, 'Quality', 90); %30x realtime
open(writerObj);

figure; title('Reposition'); pause;
subplot(4,1,1:3);
hfig = imshow(gray2rgb(handles.img0, INTENSITY_LIM));
% hfig = imagesc(handles.img0, INTENSITY_LIM); 
set(gca, {'XTick', 'YTick'}, {[],[]}); 
axis equal; axis tight; hold on; 
title('EOD rate at the head (red: higher rate)');

%plot locations
nframes = size(XC,1);
% mrColor = jet(nframes);
[mrColor, vrRateSrt, vrQuantSrt] = quantile2color(RC);

%colorbar
plotColorbar(size(handles.img0), vrRateSrt, vrQuantSrt);
EODR1 = EODR(TEOD > TLIM(1) & TEOD < TLIM(2));
RLIM = [quantile(EODR1, .001), quantile(EODR1, .999)];
htext = [];
vhChevron = [];
for iframe=1:nframes
    
    %------------------
    subplot(4,1,1:3);
    frameNum = iframe + handles.FLIM(1) - 1;
    mrImg = readFrame(handles.vidobj, frameNum);    
    mrImg(~handles.MASK) = 0;   
    mrImg = gray2rgb(mrImg, INTENSITY_LIM);
    set(hfig, 'cdata', mrImg);    

    try delete(htext); catch, end;
    htext(1) = text(10, 30, sprintf('EOD (%s): %0.1f Hz', chName, RC(iframe)), ...
        'FontSize', 12, 'Color', [1 1 1]);
    htext(2) = text(10, 75, sprintf('Time: %0.1f s', TC(iframe)), ...
        'FontSize', 12, 'Color', [1 1 1]);    
    htext(3) = text(10, 120, sprintf('Frame: %d', frameNum), ...
        'FontSize', 12, 'Color', [1 1 1]);    
    
    if mod(iframe, MOV_PLOTSTEP) == 0
        vhChevron(end+1) = plotChevron(XC(iframe, 2:3), YC(iframe, 2:3), mrColor(iframe,:), 90, .3);
        if numel(vhChevron) > MOV_PLOTLEN
            delete(vhChevron(1:end-MOV_PLOTLEN));
            vhChevron(1:end-MOV_PLOTLEN) = [];
        end
    end
    
    %------------------
    subplot(4,1,4);    
    hold off;
    plot(TEOD - TC(iframe), EODR, 'k.'); hold on;
    axis([MOV_TimeWin(1), MOV_TimeWin(2), RLIM(1), RLIM(2)]);
    plot([0 0], get(gca, 'YLim'), 'r-');
    grid on;
    xlabel('Time (sec)'); 
    ylabel(sprintf('EOD (%s) Hz', chName));  
    
    colormap jet;
%     drawnow;
    try
        writeVideo(writerObj, getframe(gcf));
    catch
        disp('Movie output cancelled by user');
        close(writerObj);
        return;
    end
end

close(writerObj);
msgbox(sprintf('File written to %s', MOV_FILEOUT));