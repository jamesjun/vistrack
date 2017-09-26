function [VEL, ACC, ANG, AVEL, XHs, YHs, TC, HANG, HAVEL, Xs, Ys] = calcVelocity(handles)

LOADSETTINGS;
TC = handles.TC;
FPS = 15;

%Smooth trajectory
Xs = filtPos(handles.XC, TRAJ_NFILT, 1);
Ys = filtPos(handles.YC, TRAJ_NFILT, 1);

% XCs = filtPos(handles.XC(:,1), TRAJ_NFILT);
% YCs = filtPos(handles.YC(:,1), TRAJ_NFILT);
% XHs = filtPos(handles.XC(:,2), TRAJ_NFILT);
% YHs = filtPos(handles.YC(:,2), TRAJ_NFILT);
% X3s = filtPos(handles.XC(:,3), TRAJ_NFILT);
% Y3s = filtPos(handles.YC(:,3), TRAJ_NFILT);
% X5s = filtPos(handles.XC(:,5), TRAJ_NFILT);
% Y5s = filtPos(handles.YC(:,5), TRAJ_NFILT);
% XTs = filtPos(handles.XC(:,6), TRAJ_NFILT);
% YTs = filtPos(handles.YC(:,6), TRAJ_NFILT);

% figure; hold on; 
% plot(TC, handles.XC(:,2), 'r');
% plot(TC, XHs, 'b');
XHs = Xs(:,2);
YHs = Ys(:,2);

OXh = Xs(:,2) - Xs(:,3); %orientation
OYh = Ys(:,2) - Ys(:,3); %orientation
OXt = Xs(:,5) - Xs(:,6); %orientation
OYt = Ys(:,5) - Ys(:,6); %orientation
VX = differentiate5(Xs(:,1), 1/FPS)';
VY = differentiate5(Ys(:,1), 1/FPS)';
VN = sqrt(VX.^2 + VY.^2);
Dir = OXh .* VX + OYh .* VY;

% tail angle calculation
vrZero = zeros(size(OXh));
vrCrossP = cross([OXh, OYh, vrZero], [OXt, OYt, vrZero]);
vrCrossP = vrCrossP(:,3);
ANG = asin(vrCrossP ./ (sqrt(OXh.^2 + OYh.^2) .* sqrt(OXt.^2 + OYt.^2)));
HANG = cart2pol(OXh(:), OYh(:)); % heading angle

VEL = filtPos(VN .* sign(Dir), TRAJ_NFILT);
% ACC = differentiate5(VEL, 1/FPS);
ACC = filtPos(differentiate5(VEL, 1/FPS), TRAJ_NFILT)';
AVEL = filtPos(differentiate5(ANG, 1/FPS), TRAJ_NFILT)';
% HAVEL = filtAng(differentiate5(unwrap(HANG), 1/FPS), TRAJ_NFILT)';
HAVEL = differentiate5(unwrap(HANG), 1/FPS);

% ACCi = interp1(TC, ACC, TEOD, 'spline');


% figure;
% AX = [];
% subplot 211; plot(TC, V, '.-'); title('signed vel'); AX(1) = gca; grid on;
% subplot 212; plot(TC, VN, '-'); title('dir'); AX(2) = gca; grid on;
% linkaxes(AX, 'x');


return;

% %Interpolate
% 
% %differentiate at 100 Hz
% 
% %calc velocity
% TCi = handles.TC(1):(1/EODR_SR):handles.TC(end);
% 
% 
% Xi = interp1(handles.TC, , TCi, 'spline', 'extrap'); 
% Yi = interp1(handles.TC, filtPos(handles.YC(:,2), TRAJ_NFILT), TCi, 'spline', 'extrap'); 
% 
% VX = differentiate5(Xi, 1/EODR_SR);
% VY = differentiate5(Yi);