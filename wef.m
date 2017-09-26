function varargout = wef(vcCommand, arg1, arg2, arg3, arg4)
% wef command

if nargin<1, vcCommand='help'; end
if nargin<2, arg1=''; end
if nargin<3, arg2=''; end
if nargin<4, arg3=''; end
if nargin<5, arg4=''; end

switch vcCommand
    case 'help'
        help_();
    case {'traj', 'trajectory'}
        traj_(arg1);
    case 'info-set'
        info_set_(arg1, arg2);
    case 'make-set' %collect trials
        make_set_(arg1);
    case 'plot-lc' %plot learning curve
        plot_lc_(arg1, arg2);
    case 'plot-probe' %plot probe
        plot_probe_(arg1, arg2, arg3);        
    case 'test'
        varargout{1} = test_(arg1, arg2, arg3, arg4);
end %switch
end %func


%--------------------------------------------------------------------------
function traj_(vcFile_prm)
S_mat = load_prm_(vcFile_prm);

mnImg = S_mat.P.I0;
mlMask = S_mat.P.MASK;
mnImg(~mlMask) = mnImg(~mlMask)/4;

figure; imshow(mnImg); hold on;
plot(S_mat.HEADXY(1,:), S_mat.HEADXY(2,:));
plot(S_mat.P.CENTERPOS(1), S_mat.P.CENTERPOS(2), 'r*');
plot(S_mat.P.CENTERPOS(1), S_mat.P.CENTERPOS(2), 'r*');
set(gcf,'Name', S_mat.EODTTL);
end %func


%--------------------------------------------------------------------------
function S_mat = load_prm_(vcFile_prm)
eval(sprintf('%s;', strrep(vcFile_prm, '.m', '')));
vcFile_mat = fullfile(vcDir, vcFile);
S_mat = load(vcFile_mat);
end %func


%--------------------------------------------------------------------------
function help_()

end %func


%--------------------------------------------------------------------------
function sync_()

end %func


%--------------------------------------------------------------------------
function plot_lc_(vcSet, vcAnimals)

% collect file names
[S_dataset, P] = load_set_(vcSet, vcAnimals);
if isempty(S_dataset), return; end
csAnimals = S_dataset.csAnimals;

mrPath_meter = pool_pathlen_(S_dataset, csAnimals, P) / 100;
mrPath_iqr = quantile(mrPath_meter, [.25,.5,.75])';

assignWorkspace_(mrPath_meter, mrPath_iqr);

figure; 
errorbar_jjj([], mrPath_iqr); xlabel('Session #'); ylabel('Dist (m)');
grid on;
set(gcf,'Name',vcSet,'Color','w');
end %func


%--------------------------------------------------------------------------
function cs = vc2cs_(vc)
cs = arrayfun(@(x)x, vc, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function mrPath = pool_pathlen_(S_dataset, csAnimals, P)
nTrialsPerSession = get_set_(P, 'nTrialsPerSession', 4);
cvrPath = cell(numel(csAnimals), 1);
for iAnimal = 1:numel(csAnimals)
%     S_dataset.vsTrial_A
    csTrial_ = getfield(S_dataset, sprintf('vsTrial_%s', csAnimals{iAnimal}));
%     csTrials_ = reshape_(getfield(S_dataset, vcName_), nTrialsPerSession);
    cvrPath{iAnimal} = reshape_(get_(csTrial_, 'pathLen_cm'), nTrialsPerSession);
end
nSessions = min(cellfun(@(x)size(x,2), cvrPath));
cvrPath = cellfun(@(x)x(:,1:nSessions), cvrPath, 'UniformOutput', 0);
mrPath = cell2mat(cvrPath);
end %func


%--------------------------------------------------------------------------
function ccPath = pool_(S_dataset, csAnimals, vcName)
ccPath = cell(numel(csAnimals), 1);
for iAnimal = 1:numel(csAnimals)
    csTrial_ = getfield(S_dataset, sprintf('vsTrial_%s', csAnimals{iAnimal}));
    ccPath{iAnimal} = cellstruct_get_(csTrial_, vcName);
end
end %func


%--------------------------------------------------------------------------
% create a dataset
function S = make_set_(vcSet)
[S_dataset, P] = load_set_(vcSet);
csAnimals = S_dataset.csAnimals;
% collect file names
% vcFile_trial_ = csFiles_trial{1};
% S_calib = calibrate_(vcFile_trial_); % correction factor for loading trial

% filter animals
S_set = struct();
for iAnimal = 1:numel(csAnimals) % go by animals and call importTrial    
    vcAnimal_ = csAnimals{iAnimal};
    vcTrial_ = sprintf('vsTrial_%s', vcAnimal_);
    vcProbe_ = sprintf('vsProbe_%s', vcAnimal_);
    [S_set.(vcTrial_), S_set.(vcProbe_)] = import_trials_(csFiles_trial, vcAnimal_);   
end

% write to file set
struct_save_(S_set, S_dataset.vcDataset);
end %func


%--------------------------------------------------------------------------
function P = mfile2struct_(vcFile_input_exclude_later_m)
% James Jun 2017 May 23
% Run a text file as .m script and result saved to a struct P
% _prm and _prb can now be called .prm and .prb files

eval(sprintf('%s;', strrep(vcFile_input_exclude_later_m, '.m', '')));
S_ws = whos(); 
csVars = {S_ws.name};
csVars = setdiff(csVars, 'vcFile_input_exclude_later_m');
P = struct();
for iVar=1:numel(csVars)
    try
        vcVar_ = csVars{iVar};
        eval(sprintf('P.(''%s'') = %s;', vcVar_, vcVar_));
    catch
        disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
% function disperr_()
% Display error message and the error stack
function disperr_(vcMsg)
% ask user to email jrclust@vidriotech.com ? for the error ticket?
dbstack('-completenames'); % display an error stack
vcErr = lasterr();
if nargin==0
    fprintf(2, '%s\n', vcErr);
else
    fprintf(2, '%s:\n\t%s\n', vcMsg, vcErr);
end
try gpuDevice(1); disp('GPU device reset'); catch, end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and added test ouput
function S_out = test_(vcFunc, cell_Input, nOutput, fVerbose)
% S_out = test_(vcFunc, {input1, input2, ...}, nOutput)

if nargin<2, cell_Input = {}; end
if nargin<3, nOutput = []; end
if nargin<4, fVerbose = ''; end

if isempty(nOutput), nOutput = 1; end
if ~iscell(cell_Input), cell_Input = {cell_Input}; end
if isempty(fVerbose), fVerbose = 1; end

delete_empty_files_(); % delete empty files

try
    switch nOutput
        case 0
            feval(vcFunc, cell_Input{:});
            S_out = [];
        case 1
            [out1] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1);
        case 2
            [out1, out2] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2);
        case 3
            [out1, out2, out3] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2, out3);
        case 4
            [out1, out2, out3, out4] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2, out3, out4);
    end %switch
    if fVerbose
        if nOutput>=1, fprintf('[%s: out1]\n', vcFunc); disp(S_out.out1); end
        if nOutput>=2, fprintf('[%s: out2]\n', vcFunc); disp(S_out.out2); end
        if nOutput>=3, fprintf('[%s: out3]\n', vcFunc); disp(S_out.out3); end
        if nOutput>=4, fprintf('[%s: out4]\n', vcFunc); disp(S_out.out4); end
    end
catch
    disperr_();
    S_out = [];
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function delete_empty_files_(vcDir)
if nargin<1, vcDir=[]; end
delete_files_(find_empty_files_(vcDir));
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function delete_files_(csFiles, fVerbose)
% Delete list of files
% delete_files_(vcFile)
% delete_files_(csFiles)
% delete_files_(csFiles, fVerbose)

if nargin<2, fVerbose = 1; end
if ischar(csFiles), csFiles = {csFiles}; end
for iFile = 1:numel(csFiles)
    try
        if exist(csFiles{iFile}, 'file')
            delete(csFiles{iFile});
            if fVerbose
                fprintf('\tdeleted %s.\n', csFiles{iFile});
            end
        end
    catch
        disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and testing
function csFiles = find_empty_files_(vcDir)
% find files with 0 bytes

if nargin==0, vcDir = []; end
if isempty(vcDir), vcDir = pwd(); end
vS_dir = dir(vcDir);
viFile = find([vS_dir.bytes] == 0 & ~[vS_dir.isdir]);
csFiles = {vS_dir(viFile).name};
csFiles = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name

S=[];
for i=1:nargin
    S = setfield(S, inputname(i), varargin{i});
end
end %func


%--------------------------------------------------------------------------
function [csFiles_full, csFiles] = dir_(vcFilter, csExcl)
% return name of files full path, exclude files
if nargin>=2
    if ischar(csExcl), csExcl = {csExcl}; end
    csExcl = union(csExcl, {'.', '..'}); 
else
    csExcl = [];
end
csFiles = dir(vcFilter);
csFiles  = {csFiles.('name')};
csFiles = setdiff(csFiles, csExcl);
[vcDir, ~, ~] = fileparts(vcFilter);
if isempty(vcDir), vcDir='.'; end
csFiles_full = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function struct_save_(S, vcFile, fVerbose)
% 7/13/17 JJJ: Version check routine
if nargin<3, fVerbose = 0; end
if fVerbose
    fprintf('Saving a struct to %s...\n', vcFile); t1=tic;
end
version_year = version('-release');
version_year = str2double(version_year(1:end-1));
if version_year >= 2017
    save(vcFile, '-struct', 'S', '-v7.3', '-nocompression'); %faster    
else
%     disp('Saving with -nocompression flag failed. Trying without compression');
    save(vcFile, '-struct', 'S', '-v7.3');
end
if fVerbose
    fprintf('\ttook %0.1fs.\n', toc(t1));
end
end %func


%--------------------------------------------------------------------------
function S_calib = calibrate_(vcFile_trial)
% determine rotation and center
S = load(vcFile_trial);
hFig = figure; 
imshow(imadjust(S.img0));
title('click (-50,0), (+50,0)cm');
set(gcf, 'Position', get(0, 'ScreenSize'));
[vrX, vrY] = ginput(2);
xy0 = [mean(vrX), mean(vrY)];
angXaxis = rad2deg(cart2pol(diff(vrX), diff(vrY))); %in rad
pixpercm = sqrt(diff(vrX)^2 + diff(vrY)^2) / 100;

hold on;
plot(vrX, vrY, 'r.');
plot(xy0(1), xy0(2), 'r.');
plot(S.xy0(1), S.xy0(2), 'go');
vcDisp = sprintf('%s, ang: %0.3f deg, pixpercm: %0.3f, x0: %0.1f, y0: %0.1f', ...
    vcFile_trial, angXaxis, pixpercm, xy0(1), xy0(2));
disp(vcDisp);
title(vcDisp);
uiwait(msgbox('press okay to continue'));
close(hFig); drawnow;

S_calib = makeStruct_(angXaxis, pixpercm, xy0);
end %func


%--------------------------------------------------------------------------
function [csTrial, csProbe] = import_trials_(csFiles, vcAnimal)
[csTrial, csProbe] = deal({});
warning off;
for iFile = 1:numel(csFiles)
    vcFile_ = csFiles{iFile};
    try        
        [fAnimal, fProbe] = checkFile_(vcFile_, vcAnimal);
        if ~fAnimal, continue; end
        if fProbe
            csProbe{end+1} = importTrial(vcFile_);
        else
            csTrial{end+1} = importTrial(vcFile_);
        end
        fprintf('%d/%d: %s\n', iFile, numel(csFiles), vcFile_);
    catch
        disperr_(vcFile_);
    end
end
end %func


%--------------------------------------------------------------------------
function [fAnimal, fProbe] = checkFile_(vcFile, vcAnimal)
% returns if the animal name matches and whether probe or not
% file name takes "%##%#'
[~, vcDataId, ~] = fileparts(vcFile);
vcDataId = strrep(vcDataId, '_Track', '');
fAnimal = upper(vcDataId(4)) == upper(vcAnimal);
fProbe = numel(vcDataId) > 5;
end %func


%--------------------------------------------------------------------------
function val = get_set_(S, vcName, def_val)
% set a value if field does not exist (empty)
if isempty(S), val = def_val; return; end
if ~isstruct(S)
    val = []; 
    fprintf(2, 'get_set_: %s be a struct\n', inputname(1));
    return;
end
val = get_(S, vcName);
if isempty(val), val = def_val; end
end %func


%--------------------------------------------------------------------------
function cs = cellstruct_get_(cS, vcName)
% return cell of info
cs = cell(size(cS));
for i=1:numel(cs)
    try
        cs{i} = cS{i}.(vcName);
    catch
        ;
    end
end 
end %func

%--------------------------------------------------------------------------
function varargout = get_(varargin)
% retrieve a field. if not exist then return empty
% [val1, val2] = get_(S, field1, field2, ...)

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end
if iscell(S)
    out1 = [];
    for i=1:numel(S)
        try
            out1(end+1) = S{i}.(varargin{2});
        catch
            ;
        end
    end
    varargout{1} = out1;
    return;
end

for i=2:nargin
    vcField = varargin{i};
    try
        varargout{i-1} = S.(vcField);
    catch
        varargout{i-1} = [];
    end
end
end %func



%--------------------------------------------------------------------------
function mr = reshape_(vr, nwin)
nbins = floor(numel(vr)/nwin);
mr = reshape(vr(1:nbins*nwin), nwin, nbins);
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function assignWorkspace_(varargin)
% Assign variables to the Workspace

for i=1:numel(varargin)
    if ~isempty(varargin{i})
        assignin('base', inputname(i), varargin{i});
        fprintf('assigned ''%s'' to workspace\n', inputname(i));
    end
end
end %func



%--------------------------------------------------------------------------
function plot_probe_(vcSet, vcAnimals, vcProbe)
[S_dataset, P] = load_probe_(vcSet, vcAnimals, vcProbe);
csAnimals = S_dataset.csAnimals;

img0 = vsTrialPool_P{1}.img0;
mrImg_P = calcVisitCount(vsTrialPool_P, img0);
mlMask = getImageMask(img0, [0 60], 'CENTRE');
figure; 
imshow(rgbmix(imadjust(img0), mrImg_P, mlMask)); 
title(sprintf('%s: %s', vcSet, sprintf('%s, ', csAnimals{:})));
grid on;
set(gcf,'Name',vcSet,'Color','w');
disp(cellstruct_get_(vsTrialPool_P, 'dataID')');
end %func


%--------------------------------------------------------------------------
function info_set_(vcSet, vcAnimals)
% vcSet: {'rand', 'randwide', 'none', 'cue', 'stable', 'shuffle'}
P = mfile2struct_('settings_wef.m');
if isempty(vcSet)
    fprintf(2, 'Specify set name:\n\t%s\n', sprintf('%s, ', P.csNames_set{:}));
    return; 
end

vcSet = lower(vcSet);
eval(sprintf('vcDir = P.vcDir_%s;', vcSet));
eval(sprintf('vcDataset = P.vcDataset_%s;', vcSet));
eval(sprintf('csAnimals = P.csAnimals_%s;', vcSet));
if ~isempty(vcAnimals), csAnimals = vc2cs_(vcAnimals); end

% Show files in the set. trial duration and day
S_dataset = load(vcDataset);
ccDataId = pool_(S_dataset, csAnimals, 'dataID');
ccvtEod = pool_(S_dataset, csAnimals, 'TEOD');
csLine = {'b.-', 'r.-', 'g.-', 'k.-'};
figure; hold on;
for iAnimal = 1:numel(ccDataId)
    csDataId_ = ccDataId{iAnimal};
    viSession_ = cellfun(@(x)str2double(x([2,3,5])), csDataId_);
    vt_dur_ = cellfun(@(x)diff(x([1,end])), ccvtEod{iAnimal});
    plot(viSession_, vt_dur_, csLine{iAnimal});
end
legend(csAnimals);
end %func


%--------------------------------------------------------------------------
function [cvtDur, cviSession, csAnimals] = pool_duration_(S_dataset, csAnimals)
if nargin<2, csAnimals = {'A', 'B', 'C', 'D'}; end

ccDataId = pool_(S_dataset, csAnimals, 'dataID');
ccvtEod = pool_(S_dataset, csAnimals, 'TEOD');
for iAnimal = 1:numel(ccDataId)
    cvtDur{iAnimal} = cellfun(@(x)diff(x([1,end])), ccvtEod{iAnimal});
    cviSession{iAnimal} = cellfun(@(x)str2double(x([2,3,5])), ccDataId{iAnimal});    
end
end %func


%--------------------------------------------------------------------------
function S_dataset = filter_duration_(S_dataset, maxDur)
[cvtDur, cviSession, csAnimals] = pool_duration_(S_dataset);
for iAnimal = 1:numel(cvtDur)
    vcField_ = sprintf('vsTrial_%s', csAnimals{iAnimal});
    csTrial_ = getfield(S_dataset, vcField_);
    S_dataset.(vcField_) = csTrial_(cvtDur{iAnimal} < maxDur);
end
end %func


%--------------------------------------------------------------------------
function S_dataset = filter_valid_(S_dataset, csAnimals)
for iAnimal = 1:numel(csAnimals)
    vcField_ = sprintf('vsTrial_%s', csAnimals{iAnimal});
    csTrial_ = getfield(S_dataset, vcField_);
    vl_ = cellfun(@(x)isfield(x, 'TEOD'), csTrial_);
    S_dataset.(vcField_) = csTrial_(vl_);
end
end %func


%--------------------------------------------------------------------------
function [S_dataset, P] = load_set_(vcSet, vcAnimals)
% vcSet: {'rand', 'randwide', 'none', 'cue', 'stable', 'shuffle'}
if nargin<2, vcAnimals = ''; end
[S_dataset, P] = deal([]);
P = mfile2struct_('settings_wef.m');

if isempty(vcSet)
    fprintf(2, 'Specify set name:\n\t%s\n', sprintf('%s, ', P.csNames_set{:}));
    return; 
end

vcSet = lower(vcSet);
eval(sprintf('vcDir = P.vcDir_%s;', vcSet)); % not needed
eval(sprintf('vcDataset = P.vcDataset_%s;', vcSet));

if ~isempty(vcAnimals)
    csAnimals = vc2cs_(vcAnimals);
else
    eval(sprintf('csAnimals = P.csAnimals_%s;', vcSet));
end

S_dataset = load(vcDataset);
S_dataset = filter_valid_(S_dataset, P.csAnimals);
S_dataset = filter_duration_(S_dataset, P.maxDur);

if vcDir(end) ~= filesep(), vcDir(end+1) = filesep(); end
csFiles_trial = dir_(sprintf('%s*_Track.mat', vcDir));

% add to the struct
S_dataset.csFiles_trial = csFiles_trial;
S_dataset.vcDataset = vcDataset;
S_dataset.csAnimals = csAnimals;

end %func


%--------------------------------------------------------------------------
function [csTrials_probe, P] = load_probe_(vcSet, vcAnimals, vcProbe)
% vcSet: {'rand', 'randwide', 'none', 'cue', 'stable', 'shuffle'}
if isempty(vcProbe), vcProbe = 'probe'; end

[csTrials_probe, P] = deal([]);
P = mfile2struct_('settings_wef.m');
if isempty(vcSet)
    fprintf(2, 'Specify set name:\n\t%s\n', sprintf('%s, ', P.csNames_set{:}));
    return; 
end
vcSet = lower(vcSet);
eval(sprintf('vcDir = P.vcDir_%s;', vcSet)); % not needed
eval(sprintf('vcDataset = P.vcDataset_%s;', vcSet));

if ~isempty(vcAnimals)
    csAnimals = vc2cs_(vcAnimals);
else
    eval(sprintf('csAnimals = P.csAnimals_%s;', vcSet));
end

% collect trials directly
vcDir_relearn
S_dataset = load(vcDataset);
S_dataset = filter_valid_(S_dataset, P.csAnimals);
S_dataset = filter_duration_(S_dataset, P.maxDur);

if vcDir(end) ~= filesep(), vcDir(end+1) = filesep(); end
csFiles_trial = dir_(sprintf('%s*_Track.mat', vcDir));

% add to the struct
S_dataset.csFiles_trial = csFiles_trial;
S_dataset.vcDataset = vcDataset;
S_dataset.csAnimals = csAnimals;

end %func

