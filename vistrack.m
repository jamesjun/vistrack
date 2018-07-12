function varargout = vistrack(varargin)
vcCmd = 'help';
if nargin>=1, vcCmd = varargin{1}; else vcCmd = 'help'; end
if nargin>=2, vcArg1 = varargin{2}; else vcArg1 = ''; end
if nargin>=3, vcArg2 = varargin{3}; else vcArg2 = ''; end
if nargin>=4, vcArg3 = varargin{4}; else vcArg3 = ''; end
if nargin>=5, vcArg4 = varargin{5}; else vcArg4 = ''; end
if nargin>=6, vcArg5 = varargin{6}; else vcArg5 = ''; end

% Command interpreter
fReturn = 1;
switch lower(vcCmd)
    case 'commit',  commit_(); 
    case 'help', help_(vcArg1); 
    case 'issue', issue_(vcArg1); 
    case 'wiki', wiki_(vcArg1); 
    case 'version'
        if nargout>0
            [varargout{1}, varargout{2}] = version_();
        else
            version_();
        end
        return;
    case 'gui', GUI(); 
    case 'edit', edit_(vcArg1); 
    case 'unit-test', unit_test_(vcArg1);    
    case 'update', update_(vcArg1);
    
    case 'measure_trials', [varargout{1}, varargout{2}] = measure_trials_(vcArg1, vcArg2);
        
    otherwise, help_();
end %switch
if fReturn, return; end

end %func


%--------------------------------------------------------------------------
function commit_(vcDir_target)
if nargin<1, vcDir_target = 'D:\Dropbox\Git\vistrack\'; end

% Delete previous files
S_warning = warning();
warning('off');
delete_empty_files_();
delete([vcDir_target, '*']);
warning(S_warning);

% Copy files
csFiles_upload = {'*.m', 'GUI.fig', 'change_log.txt', 'readme.txt', 'example.trialset'};
copyfile_(csFiles_upload, vcDir_target, '.');

edit_('change_log.txt');
end %func


%--------------------------------------------------------------------------
function delete_empty_files_(vcDir)
if nargin<1, vcDir=[]; end
delete_files_(find_empty_files_(vcDir));
end %func


%--------------------------------------------------------------------------
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
function delete_files_(csFiles)
for iFile = 1:numel(csFiles)
    try
        if exist(csFiles{iFile}, 'file')
            delete(csFiles{iFile});
            fprintf('\tdeleted %s.\n', csFiles{iFile});
        end
    catch
        disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
% 9/29/17 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate] = version_()
if nargin<1, vcFile_prm = ''; end
vcVer = 'v0.1.6';
vcDate = '7/12/2018';
if nargout==0
    fprintf('%s (%s) installed\n', vcVer, vcDate);
end
end %func


%--------------------------------------------------------------------------
function csHelp = help_(vcCommand)
if nargin<1, vcCommand = ''; end
if ~isempty(vcCommand), wiki_(vcCommand); return; end

csHelp = {...
    ''; 
    'Usage: vistrack command arg1 arg2 ...';
    '';
    '# Documentation';
    '  vistrack help';
    '    Display a help menu';       
    '  vistrack version';
    '    Display the version number and the updated date';                     
    '  vistrack wiki';
    '    Open vistrack Wiki on GitHub';     
    '  vistrack issue';
    '    Post an issue at GitHub (log-in with your GitHub account)';        
    '';
    '# Main commands';
    '  vistrack edit (mysettings.prm)';
    '    Edit .vistrack file currently working on'; 
    '  vistrack setprm myparam.prm';
    '    Select a .prm file to use'; 
    '  vistrack clear';
    '    Clear cache';
    '  vistrack clear myparam.prm';
    '    Delete previous results (files: _jrc.mat, _spkwav.jrc, _spkraw.jrc, _spkfet.jrc)';        
    '';
    '# Batch process';
    '  vistrack dir myparam.prm'; 
    '    List all recording files to be clustered together (csFile_merge)';
    '';
    '# Deployment';
    '  vistrack update';
    '    Update from Github'; 
    '  vistrack update version';
    '    Download specific version from Github'; 
    '  vistrack commit';
    '    Commit vistrack code to Github';
    '  vistrack unit-test';
    '    Run a suite of unit teste.';       
    '';
};
if nargout==0, disp_cs_(csHelp); end
end %func


%--------------------------------------------------------------------------
function disp_cs_(cs)
% display cell string
cellfun(@(s)fprintf('%s\n',s), cs);
end %func


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function issue_(vcMode)
% issue_
% issue_ post
if nargin<1, vcMode = 'search'; end
switch lower(vcMode)
    case 'post', web_('https://github.com/jamesjun/vistrack/issues/new')
    otherwise, web_('https://github.com/jamesjun/vistrack/issues')
end %switch
end %func


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function wiki_(vcPage)
if nargin<1, vcPage = ''; end
if isempty(vcPage)
    web_('https://github.com/jamesjun/vistrack/wiki'); 
else
    web_(['https://github.com/jamesjun/vistrack/wiki/', vcPage]); 
end
end %func


%--------------------------------------------------------------------------
function web_(vcPage)
if isempty(vcPage), return; end
if ~ischar(vcPage), return; end
try
    % use system browser
    if ispc()
        system(['explorer ', vcPage]);
    elseif ismac()
        system(['open ', vcPage]);
    elseif isunix()
        system(['gnome-open ', vcPage]);
    else
        web(vcPage);
    end
catch
    web(vcPage); % use matlab default web browser
end
end %func


%--------------------------------------------------------------------------
% 10/8/17 JJJ: Created
% 3/20/18 JJJ: captures edit failure (when running "matlab -nodesktop")
function edit_(vcFile)
% vcFile0 = vcFile;
if isempty(vcFile), vcFile = mfilename(); end
if ~exist_file_(vcFile)
    fprintf(2, 'File does not exist: %s\n', vcFile);
    return; 
end
fprintf('Editing %s\n', vcFile);
try edit(vcFile); catch, end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile)
if nargin<2, fVerbose = 0; end
if isempty(vcFile)
    flag = 0; 
else
    flag = ~isempty(dir(vcFile));
end
if fVerbose && ~flag
    fprintf(2, 'File does not exist: %s\n', vcFile);
end
end %func


%--------------------------------------------------------------------------
function nFailed = unit_test_(vcArg1, vcArg2, vcArg3)
% 2017/2/24. James Jun. built-in unit test suite (request from Karel Svoboda)
% run unit test
%[Usage]
% unit_test()
%   run all
% unit_test(iTest)
%   run specific test again and show profile
% unit_test('show')
%   run specific test again and show profile
% @TODO: test using multiple datasets and parameters.
global fDebug_ui;

if nargin<1, vcArg1 = ''; end
if nargin<2, vcArg2 = ''; end
if nargin<3, vcArg3 = ''; end

cd(fileparts(mfilename('fullpath'))); % move to jrclust folder
% if ~exist_file_('sample.bin'), jrc3('download', 'sample'); end

nFailed = 0;
profile('clear'); %reset profile stats
csCmd = {...  
    'close all; clear all;', ... %start from blank
    'vistrack', ...
    'vistrack help', ...
    'vistrack version', ...
    'vistrack wiki', ...
    'vistrack issue', ...
    }; %last one should be the manual test

if ~isempty(vcArg1)
    switch lower(vcArg1)
        case {'show', 'info', 'list', 'help'}
            arrayfun(@(i)fprintf('%d: %s\n', i, csCmd{i}), 1:numel(csCmd)); 
            return;
        case {'manual', 'ui', 'ui-manual'}
            iTest = numel(csCmd); % + [-1,0];
        case {'traces', 'ui-traces'}
            iTest = numel(csCmd)-2; % second last
        otherwise
            iTest = str2num(vcArg1);
    end        
    fprintf('Running test %s: %s\n', vcArg1, csCmd{iTest});
    csCmd = csCmd(iTest);
end

vlPass = false(size(csCmd));
[csError, cS_prof] = deal(cell(size(csCmd)));
vrRunTime = zeros(size(csCmd));
for iCmd = 1:numel(csCmd)   
    eval('close all; fprintf(''\n\n'');'); %clear memory
    fprintf('Test %d/%d: %s\n', iCmd, numel(csCmd), csCmd{iCmd}); 
    t1 = tic;        
    profile('on');
    fDebug_ui = 1;
%     set0_(fDebug_ui);
    try
        if any(csCmd{iCmd} == '(' | csCmd{iCmd} == ';') %it's a function        
            evalin('base', csCmd{iCmd}); %run profiler 
        else % captured by profile
            csCmd1 = strsplit(csCmd{iCmd}, ' ');
            feval(csCmd1{:});
        end
        vlPass(iCmd) = 1; %passed test        
    catch
        csError{iCmd} = lasterr();
        fprintf(2, '\tTest %d/%d failed\n', iCmd, numel(csCmd));
    end
    vrRunTime(iCmd) = toc(t1);
    cS_prof{iCmd} = profile('info');    
end
nFailed = sum(~vlPass);

fprintf('Unit test summary: %d/%d failed.\n', sum(~vlPass), numel(vlPass));
for iCmd = 1:numel(csCmd)
    if vlPass(iCmd)
        fprintf('\tTest %d/%d (''%s'') took %0.1fs.\n', iCmd, numel(csCmd), csCmd{iCmd}, vrRunTime(iCmd));
    else
        fprintf(2, '\tTest %d/%d (''%s'') failed:%s\n', iCmd, numel(csCmd), csCmd{iCmd}, csError{iCmd});
    end        
end

if numel(cS_prof)>1
    assignWorkspace_(cS_prof);
    disp('To view profile, run: profview(0, cS_prof{iTest});');
else
    profview(0, cS_prof{1});
end
fDebug_ui = [];
% set0_(fDebug_ui);
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Output message is added
% 8/2/17 JJJ: Test and documentation
function vcMsg = assignWorkspace_(varargin)
% Assign variables to the Workspace
vcMsg = {};
for i=1:numel(varargin)
    if ~isempty(varargin{i})
        assignin('base', inputname(i), varargin{i});
        vcMsg{end+1} = sprintf('assigned ''%s'' to workspace\n', inputname(i));        
    end
end
vcMsg = cell2mat(vcMsg);
if nargout==0, fprintf(vcMsg); end
end %func


%--------------------------------------------------------------------------
function update_(vcVersion)
if ~exist_dir_('.git')
    fprintf(2, 'Not a git repository. run "git clone https://github.com/jamesjun/vistrack"\n');
    return; 
end
if nargin<1, vcVersion = ''; end
delete_('settings_vistrack.m');
delete_('example.trialset');
repoURL = 'https://github.com/jamesjun/vistrack';
try
    if isempty(vcVersion)
        code = system('git pull');
    else
        code = system(sprintf('git reset --hard "%s"', vcVersion));
    end
catch
	code = -1;
end
if code ~= 0
    fprintf(2, 'Not a git repository. Please run the following command to clone from GitHub.\n');    
    fprintf(2, '\tRun system(''git clone %s.git''\n', repoURL);
    fprintf(2, '\tor install git from https://git-scm.com/downloads\n');  
else
    edit('change_log.txt');
end
end %func


%--------------------------------------------------------------------------
% 11/5/17 JJJ: added vcDir_from
% 9/26/17 JJJ: multiple targeting copy file. Tested
function copyfile_(csFiles, vcDir_dest, vcDir_from)
% copyfile_(vcFile, vcDir_dest)
% copyfile_(csFiles, vcDir_dest)
% copyfile_(csFiles, csDir_dest)

if nargin<3, vcDir_from = ''; end
% Recursion if cell is used
if iscell(vcDir_dest)
    csDir_dest = vcDir_dest;
    for iDir = 1:numel(csDir_dest)
        try
            copyfile_(csFiles, csDir_dest{iDir});
        catch
            disperr_();
        end
    end
    return;
end

if ischar(csFiles), csFiles = {csFiles}; end
for iFile=1:numel(csFiles)
    vcPath_from_ = csFiles{iFile};   
    if ~isempty(vcDir_from), vcPath_from_ = fullfile(vcDir_from, vcPath_from_); end
    if exist_dir_(vcPath_from_)
        [vcPath_,~,~] = fileparts(vcPath_from_);
        vcPath_from_ =  sprintf('%s%s*', vcPath_, filesep());
        vcPath_to_ = sprintf('%s%s%s', vcDir_dest, filesep(), dir_filesep_(csFiles{iFile}));
        mkdir_(vcPath_to_);    
%         disp([vcPath_from_, '; ', vcPath_to_]);
    else
        vcPath_to_ = vcDir_dest;
        fCreatedDir_ = mkdir_(vcPath_to_); 
        if fCreatedDir_          
            disp(['Created a folder ', vcPath_to_]);
        end
    end    
    try   
        vcEval1 = sprintf('copyfile ''%s'' ''%s'' f;', vcPath_from_, vcPath_to_);
        eval(vcEval1);
        fprintf('\tCopied ''%s'' to ''%s''\n', vcPath_from_, vcPath_to_);
    catch
        fprintf(2, '\tFailed to copy ''%s''\n', vcPath_from_);
    end
end
end %func


%--------------------------------------------------------------------------
% 11/5/17 JJJ: Created
function flag = exist_dir_(vcDir)
if isempty(vcDir)
    flag = 0;
else
    flag = exist(vcDir, 'dir') == 7;
end
end %func


%--------------------------------------------------------------------------
function fCreatedDir = mkdir_(vcDir)
% make only if it doesn't exist. provide full path for dir
fCreatedDir = exist_dir_(vcDir);
if ~fCreatedDir
    try
        mkdir(vcDir); 
    catch
        fCreatedDir = 0;
    end
end
end %func


%--------------------------------------------------------------------------
% function xyLED = fixLedPos_(vidobj, FLIM1)
% if nargin<2, FLIM1 = []; end
% square_len = 50;
% 
% if isempty(FLIM1), FLIM1 = [1, 300]; end % show first 300 frames
% FLIM1(2) = min(vidobj.NumberOfFrames, FLIM1(2));
% try
%     hMsg = msgbox('Loading video...');
%     trImg = read(vidobj, FLIM1);
%     close_(hMsg);
% catch
%     disperr_();
%     xyLed = [];
%     return;
% end
% trImg = squeeze(trImg(:,:,1,:));   % extract red color
% % img_mean = uint8(mean(single(trImg),3));
% % img_sd = (std(single(trImg),1,3));
% img_pp = (max(trImg,[],3) - min(trImg,[],3));
% [~,imax_pp] = max(img_pp(:));
% [yLed, xLed] = ind2sub(size(img_pp), imax_pp);
% 
% hFig = create_figure_('fig_led', [0 0 1 1]); hold on;
% % subplot 121; imshow(img_mean); axis equal;
% imagesc(img_pp); axis tight; axis equal; colormap gray;
% hold on; plot(xLed, yLed, 'ro');
% % h_point = imrect(gca, [xLed, yLed, 0, 0] + square_len*[-.5,-.5,1,1])
% % h_point.setColor('r');
% 
% % xyLED = ginput(1);
% end %func


%--------------------------------------------------------------------------
% 17/12/5 JJJ: Error info is saved
% Display error message and the error stack
function disperr_(vcMsg, hErr)
% disperr_(vcMsg): error message for user
% disperr_(vcMsg, hErr): hErr: MException class
% disperr_(vcMsg, vcErr): vcErr: error string
try
    dbstack('-completenames'); % display an error stack
    if nargin<1, vcMsg = ''; end
    if nargin<2, hErr = lasterror('reset');  end
    if ischar(hErr) % properly formatted error
        vcErr = hErr;
    else
%         save_err_(hErr, vcMsg); % save hErr object?   
        vcErr = hErr.message;        
    end
catch
    vcErr = '';
end
if nargin==0
    fprintf(2, '%s\n', vcErr);
elseif ~isempty(vcErr)
    fprintf(2, '%s:\n\t%s\n', vcMsg, vcErr);
else
    fprintf(2, '%s:\n', vcMsg);
end
% try gpuDevice(1); disp('GPU device reset'); catch, end
end %func


%--------------------------------------------------------------------------
function hFig = create_figure_(vcTag, vrPos, vcName, fToolbar, fMenubar)
if nargin<2, vrPos = []; end
if nargin<3, vcName = ''; end
if nargin<4, fToolbar = 0; end
if nargin<5, fMenubar = 0; end
if isempty(vcTag)
    hFig = figure();
elseif ischar(vcTag)
    hFig = figure_new_(vcTag); 
else
    hFig = vcTag;
end
set(hFig, 'Name', vcName, 'NumberTitle', 'off', 'Color', 'w');
clf(hFig);
set(hFig, 'UserData', []); %empty out the user data
if ~fToolbar
    set(hFig, 'ToolBar', 'none'); 
else
    set(hFig, 'ToolBar', 'figure'); 
end
if ~fMenubar
    set(hFig, 'MenuBar', 'none'); 
else
    set(hFig, 'MenuBar', 'figure'); 
end

if ~isempty(vrPos), resize_figure_(hFig, vrPos); end
end %func


%--------------------------------------------------------------------------
function close_(hMsg)
try close(hMsg); catch; end
end %func


%--------------------------------------------------------------------------
function hFig = figure_new_(vcTag)
%remove prev tag duplication
delete_multi_(findobj('Tag', vcTag, 'Type', 'Figure')); 

hFig = figure('Tag', vcTag);
end %func


%--------------------------------------------------------------------------
function hFig = resize_figure_(hFig, posvec0, fRefocus)
if nargin<3, fRefocus = 1; end
height_taskbar = 40;

pos0 = get(groot, 'ScreenSize'); 
width = pos0(3); 
height = pos0(4) - height_taskbar;
% width = width;
% height = height - 132; %width offset
% width = width - 32;
posvec = [0 0 0 0];
posvec(1) = max(round(posvec0(1)*width),1);
posvec(2) = max(round(posvec0(2)*height),1) + height_taskbar;
posvec(3) = min(round(posvec0(3)*width), width);
posvec(4) = min(round(posvec0(4)*height), height);
% drawnow;
if isempty(hFig)
    hFig = figure; %create a figure
else
    hFig = figure(hFig);
end
drawnow;
set(hFig, 'OuterPosition', posvec, 'Color', 'w', 'NumberTitle', 'off');
end %func


%--------------------------------------------------------------------------
function delete_multi_(varargin)
% provide cell or multiple arguments
for i=1:nargin
    try
        vr1 = varargin{i};
        if numel(vr1)==1
            delete(varargin{i}); 
        elseif iscell(vr1)
            for i1=1:numel(vr1)
                try
                    delete(vr1{i1});
                catch
                end
            end
        else
            for i1=1:numel(vr1)
                try
                    delete(vr1(i1));
                catch
                end
            end
        end
    catch
    end
end
end %func


%--------------------------------------------------------------------------
function delete_(varargin)
for i=1:narging
    try
        delete(varargin{i});
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
function [cvrPathLen, cvrDuration] = measure_trials_(csFiles, csAnimals)

% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
pixpercm = 7.238; % run S141106_LearningCurve_Control.m first cell
angXaxis = -0.946;

for iAnimal = 1:numel(csAnimals)
    vcAnimal = csAnimals{iAnimal};
    [~, csFiles_, ~] = cellfun(@(x)fileparts(x), csFiles, 'UniformOutput', 0);
    vi_ = cellfun(@(x)x(4) == vcAnimal(1), csFiles_) & cellfun(@(x)lower(x(6)) ~= 'p', csFiles_);
    csFiles_ = csFiles(vi_);
    [vrPathLen_, vrDuration_] = deal(nan(size(csFiles_)));
    fprintf('\nLoading animal %s\n\t', vcAnimal);
    for iTrial = 1:numel(csFiles_)
        try
%             S_ = importTrial(csFiles_{iTrial}, angXaxis, pixpercm);
            S_ = load(csFiles_{iTrial});
            vrPathLen_(iTrial) = trial_pathlen_(S_, pixpercm, angXaxis);
            vrDuration_(iTrial) = diff(S_.TC([1,end]));
            fprintf('.');
        catch
            disp(lasterr());
        end
    end %for
    [cvrPathLen{iTrial}, cvrDuration{iTrial}] = deal(vrPathLen_, vrDuration_);
end %for
end %func


%--------------------------------------------------------------------------
function [pathLen_cm, XH, YH, TC1] = trial_pathlen_(S_trial, pixpercm, angXaxis)
[TC, XHc, YHc] = deal(S_trial.TC, S_trial.XC(:,2), S_trial.YC(:,2));
TC1 = linspace(TC(1), TC(end), numel(TC)*10);
XH = interp1(TC, XHc, TC1, 'spline');
YH = interp1(TC, YHc, TC1, 'spline');
pathLen = sum(sqrt(diff(XH).^2 + diff(YH).^2));

xyStart = trial_xyStart_(S_trial, pixpercm, angXaxis);

pathLen = pathLen + sqrt((XH(1) - xyStart(1)).^2 + (YH(1) - xyStart(2)).^2);
pathLen_cm = pathLen / pixpercm;

end %func


%--------------------------------------------------------------------------
function [xyStart, xyFood] = trial_xyStart_(S_trial, pixpercm, angXaxis)
[~, dataID, ~] = fileparts(S_trial.vidFname);
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
iAnimal = fishID - 'A' + 1;
rotMat = rotz(angXaxis);    rotMat = rotMat(1:2, 1:2);
xyStart = (xyStart * rotMat) .* [1, -1] * pixpercm + S_trial.xy0; %convert to pixel unit
xyFood = (xyFood * rotMat) .* [1, -1] * pixpercm + S_trial.xy0; %convert to pixel unit
end %func