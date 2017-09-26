%GUI_EODMOVIE

nPlot = 6;

%% Calculate timeseries
LOADSETTINGS;

%show the plot and set the range to play movie
%create an imline to doubleclick.
[EODR, TEOD, chName] = getSpike2Chan(handles.ADC, ADC_CH_EODR);
AMPL = getSpike2Chan(handles.ADC, ADC_CH_AMPL);
EODA = smoothFilter(differentiate5(EODR, .01), 5);
TLIM = handles.TC([1 end]);
IDXLIM = [];
IDXLIM(1) = find(TEOD > TLIM(1), 1, 'first');
IDXLIM(2) = find(TEOD < TLIM(2), 1, 'last');
IDXRNG = IDXLIM(1):IDXLIM(2);
TEOD = TEOD(IDXRNG);
EODR = EODR(IDXRNG);
EODA = EODA(IDXRNG);
[vrESAC, viESAC] = findESAC(EODA);
vtESAC = TEOD(viESAC);

% Kinematics calculation
[VEL, ACC, ANG, AVEL, XHs, YHs, TC] = calcVelocity(handles);

%% Plot
hFig = figure;
AX = [];
subplot(nPlot,2,2); plot(TEOD, EODR); AX(1) = gca; ylabel('EOD Rate (Hz)'); 
subplot(nPlot,2,4); plot(TEOD, EODA); AX(2) = gca; ylabel('EODA (Hz/s)'); 
subplot(nPlot,2,6); plot(TC, VEL, '-'); AX(3) = gca; ylabel('Vel. (pix/s)'); 
subplot(nPlot,2,8); plot(TC, ACC, '-'); AX(4) = gca; ylabel('Acc. (pix/s^2)'); 
subplot(nPlot,2,10); plot(TC, ANG, '-'); AX(5) = gca; ylabel('Ang. (rad)'); 
subplot(nPlot,2,12); plot(TC, AVEL, '-'); AX(6) = gca; ylabel('A.Vel. (rad/s)'); 
linkaxes(AX, 'x');
xlabel('Time (s)'); 
title('Set time range and double-click');
set(hFig, 'Name', handles.vidFname);
for i=1:nPlot
     set(AX(i), 'XLim', TLIM); 
     hold(AX(i), 'on');
     grid(AX(i), 'on');
end;
set(AX(1), 'YLim', [min(EODR), max(EODR)]);
set(AX(2), 'YLim', [min(EODA), max(EODA)]);
set(AX(3), 'YLim', [min(VEL), max(VEL)]);
set(AX(4), 'YLim', [min(ACC), max(ACC)]);
set(AX(5), 'YLim', pi/2*[-1 1]);
set(AX(6), 'YLim', [min(AVEL), max(AVEL)]);

subplot(nPlot,2,[1,3,5,7,9,11]);
hImg = imshow(imadjust(handles.MOV(:,:,round(end/2))));
AX(end+1) = gca;

%% plot the ESACCADES
%top 5% marking

threshESAC = quantile(vrESAC, .95);
viSel = find(vrESAC >= threshESAC);

plot(AX(2), get(AX(2), 'XLim'), threshESAC * [1 1], 'r-');
plot(AX(2), vtESAC(viSel), vrESAC(viSel), 'r.');

%% create a bar
chLine = {};
xpos = mean(get(AX(1), 'XLim'));
for i=1:nPlot
    chLine{i} = imline(AX(i), xpos * [1 1], get(AX(i),'YLim'));
    setColor(chLine{i}, [1 0 0]);
    setPositionConstraintFcn(chLine{i}, makeConstrainToRectFcn('imline',get(AX(i),'XLim'),get(AX(i),'YLim')));    
end
for i=1:nPlot    
    addNewPositionCallback(chLine{i}, @(pos) showMovieFrame(pos, handles, hImg, hFig, chLine, AX(end)));    
end

% show movie frame
% figure(hfig); hold on;
% subplot 121; hold on; %set(gca, 'XLim', TLIM);
% subplot 122; imshow(handles.MOV(:,:,IDXLIM(1)));

% implay(handles.MOV(:,:,IDXLIM(1):IDXLIM(2)));
% mrImg = handles.MOV(:,:,iMov);
% figure; imshow(mrImg); 
% title(sprintf('T_{ADC} = %0.3f s | F_{MOV} = %d', tsel, iMov));