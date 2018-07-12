function S = importTrial(vcFile, pixpercm, angXaxis)
if nargin<2, pixpercm = 1053.28/(sqrt(2)*100); end
if nargin<3, angXaxis = -1.1590; end %degrees

handles = load(vcFile);

eval('settings');
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
[VELc, ACCc, ANGc, AVELc, XHc, YHc, TC, HANGc, HAVELc, mrXc, mrYc] = calcVelocity(handles);

% Interpolate at the TEOD time
vxESAC = interp1(TC, filtPos(handles.XC(:,2), TRAJ_NFILT), TEOD(viESAC), 'spline');
vyESAC = interp1(TC, filtPos(handles.YC(:,2), TRAJ_NFILT), TEOD(viESAC), 'spline');
VEL = interp1(TC, VELc, TEOD, 'spline');
ACC = interp1(TC, ACCc, TEOD, 'spline');
ANG = interp1(TC, ANGc, TEOD, 'spline');
HANG = interp1(TC, unwrap(HANGc), TEOD, 'spline');
AVEL = interp1(TC, AVELc, TEOD, 'spline');
HAVEL = interp1(TC, HAVELc, TEOD, 'spline');
XH = interp1(TC, XHc, TEOD, 'spline');
YH = interp1(TC, YHc, TEOD, 'spline');
mrX = interp1mat(TC, mrXc, TEOD, 'spline');
mrY = interp1mat(TC, mrYc, TEOD, 'spline');
img0 = handles.img0;
xy0 = handles.xy0;
[~, dataID, ~] = fileparts(handles.vidFname);
fishID = dataID(4);
switch fishID
    case 'A'
        xyStart = [55, 50]; xyFood = [0 -10]; angRot = 0;
    case 'B'
        xyStart = [50, -55]; xyFood = [-10 0]; angRot = 90;
    case 'C'
        xyStart = [-55, -50]; xyFood = [0 10]; angRot = 180;
    case 'D'
        xyStart = [-50, 55]; xyFood = [10 0]; angRot = 270;
end
HANG = wrap(HANG - deg2rad(angRot)); %heading angle toward food
iAnimal = fishID - 'A' + 1;
rotMat = rotz(angXaxis);    rotMat = rotMat(1:2, 1:2);
xyStart = (xyStart * rotMat) .* [1, -1] * pixpercm + xy0; %convert to pixel unit
xyFood = (xyFood * rotMat) .* [1, -1] * pixpercm + xy0; %convert to pixel unit

%calc pathlen
pathLen = sum(sqrt(diff(XH).^2 + diff(YH).^2));
pathLen = pathLen + sqrt((XH(1) - xyStart(1)).^2 + (YH(1) - xyStart(2)).^2);
pathLen_cm = pathLen / pixpercm;

%Zone specification. incorrect due to XH not being rotated
Rfood = sqrt((XH - xyFood(1)).^2 + (YH - xyFood(2)).^2) / pixpercm;
Rcentre = sqrt((XH - xy0(1)).^2 + (YH - xy0(2)).^2) / pixpercm;
vlZone = (Rfood > 5) & (Rcentre <= 60);

%rotate XH to A's position for trajectory pooling
[vrX, vrY] = rotXY(XH, YH, xy0, angRot, 1);
[mrX, mrY] = rotXY(mrX, mrY, xy0, angRot, 1);

%quantile of EODR
EODRq = calcQuantileShift(EODR);

%quantile of EODR
EODRz = zscore(EODR);

duration = diff(TEOD([1 end]));

% EOD timestamps
try
    eval(sprintf('vtEOD = handles.ADCTS.%s_Ch10.times;', dataID));
catch
    csFields = fieldnames(handles.ADCTS);
    vi = cellstrFind(csFields, '_Ch10');
    if ~isempty(vi)
        vtEOD = getfield(handles.ADCTS, csFields{vi(1)});
        vtEOD = vtEOD.times;
    else
        error('no Ch10 TEOD found');
    end
end
viKeep = vtEOD > TEOD(1) & vtEOD < TEOD(end);
vtEOD = vtEOD(viKeep);

S = makeStruct(TEOD, EODR, EODA, EODRq, EODRz, vtEOD, ...
        XH, YH, VEL, ACC, ANG, AVEL, HANG, HAVEL, ...
        viESAC, vrESAC, vtESAC, vxESAC, vyESAC, img0, ...
        dataID, pathLen_cm, xyFood, xy0, xyStart, ...
        vlZone, Rfood, Rcentre, vrX, vrY, mrX, mrY, iAnimal, duration);
end

