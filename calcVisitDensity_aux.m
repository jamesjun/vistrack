function [VISITCNT, TIMECNT] = calcVisitDensity_aux(handles)

TRAJ_NFILT = 3;
SF = .25;

I0 = imresize(handles.img0, SF);
XC = handles.XC(:,2) * SF;
YC = handles.YC(:,2) * SF;

options = 'R = 3;';
tic
[VISITCNT, TIMECNT] = calcVisitDensity(I0, handles.TC, XC, YC, TRAJ_NFILT, options);
toc