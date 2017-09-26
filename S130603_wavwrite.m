%% wave file generation
%import timestamps and cut and embed. 1k sampling rate
S = load('C:\SkyDrive\Ms_jove\matlab\sample_TS.mat');
TLIM = [2 12];
Fs = 8000;

csFieldnames = fieldnames(S);
S = getfield(S, csFieldnames{1});
    
Teod = S.times;
Teod1 = Teod(Teod > TLIM(1) & Teod < TLIM(2));
viEod1 = round((Teod1 - TLIM(1)) * Fs);
tdur = diff(TLIM);
ns = round(tdur * Fs);

mlBinvec = zeros(ns,1);
mlBinvec(viEod1) = 1;

wavwrite(mlBinvec, Fs, 'output.wav');