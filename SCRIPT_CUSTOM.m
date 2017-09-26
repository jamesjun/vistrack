%GUI_EODMOVIE
LOADSETTINGS;
nPlot = 6;

%show the plot and set the range to play movie
%create an imline to doubleclick.
[EODR, TEOD, chName] = getSpike2Chan(handles.ADC, ADC_CH_EODR);
AMPL = getSpike2Chan(handles.ADC, ADC_CH_AMPL);
EODA = smoothFilter(differentiate5(EODR, .01), 5);
% [VEL, ACC, TC, VELi, ACCi] = calcVelocity(handles, TEOD);
[VEL, ACC, ANG, AVEL, TC] = calcVelocity(handles);

TLIM = handles.TC([1 end]);
IDXLIM = [];
IDXLIM(1) = find(TEOD > TLIM(1), 1, 'first');
IDXLIM(2) = find(TEOD < TLIM(2), 1, 'last');
IDXRNG = IDXLIM(1):IDXLIM(2);
TEOD = TEOD(IDXRNG);
EODR = EODR(IDXRNG);
EODA = EODA(IDXRNG);

