% Time sync parameters
%TLIM0 = [10, 390]; %First and last LED timing in ADC (optional)
ADC_CH_TCAM = 5; % ADC channel for TCAM (optional)
ADC_CH_TEOD = 10; % ADC channel for EOD Ts (optional)
ADC_CH_TEXT = 30; % ADC channel for EOD Ts (optional)
ADC_CH_EODR = 702; % ADC channel for EOD Ts (optional)
ADC_CH_AMPL = 703; %amplitude channel
ADC_CH_ESAC = 401; %ESAC chan
% FPS = 15; %frames per sec (optional override)
%TLIM = [317, 347]; %time range to track, [317, 634] (optional override)
%TLIM = [317, 347]; %time range to track, [317, 634] (optional override)
xyLED = [1554, 563]; % for 2013

% Image processing parameters
IM_THRESH = 20; % Initial threshold
winlen = 128*2.5; %Tracking window length in pixels
IM_MINCONTRAST = .02; %Minimum contrast difference threshold
ThreshLim = [0 30]; % Minimum and maximum adaptive intensity threshold
BW_SMOOTH = 5; % Edge smoothing size
TRACK_SHOW = 0; %show tracking result, 50% faster

% Plotting parameters
TRAJ_NFILT = 3; %Smoothing window
TRAJ_STEP = 8; %number of frames to skip
EODR_SR = 15; %Sampling rate for EOD Rate display
ADC_CH_PLOT = 702; %spike2 channel to plot
REPLAY_STEP = 4;

% Movie output parameters
INTENSITY_LIM = [40, 120]; %intensity display range
MOV_TimeWin = [-2, 2]; %time bar display range
MOV_FILEOUT = 'output.mp4'; %video output file name

% Sound output parameters
WAV_Fs = 8000;
WAV_FILEOUT = 'output.wav';