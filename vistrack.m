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
    case 'summary', varargout{1} = summary_(vcArg1);
    case 'export', export_(vcArg1);
    case 'videoreader', varargout{1} = VideoReader_(vcArg1);
    case 'dependencies', disp_dependencies_();
    case 'download-sample', download_sample_();
        
    case 'trialset-list', trialset_list_(vcArg1);    
    case 'trialset-learningcurve', trialset_learningcurve_(vcArg1);
    case 'trialset-barplots', trialset_barplots_(vcArg1);
        
    otherwise, help_();
end %switch
if fReturn, return; end

end %func


%--------------------------------------------------------------------------
function csMsg = summary_(handles)
t_dur = diff(handles.TC([1,end]));
[mrXYh_cam, vrT_cam] = get_traj_(handles);
[pathLen_cm, XH, YH, TC1] = trial_pathlen_(handles);
% [~, dataID, ~] = fileparts(handles.vidFname);
handles.vcVer = get_set_(handles, 'vcVer', 'pre v0.1.7');
handles.vcVer_date = get_set_(handles, 'vcVer_date', 'pre 7/20/2018');

vcFile_trial = get_(handles, 'editResultFile', 'String');
[vcFile_trial, S_dir] = fullpath_(vcFile_trial);
dataID = strrep(S_dir.name, '_Track.mat', '');

nDaysAgo = floor(now() - get_(S_dir, 'datenum'));

csMsg = {...
    sprintf('DataID: %s', dataID); 
    sprintf('  duration: %0.3f sec', t_dur); 
    sprintf('  path-length: %0.3f m', pathLen_cm/100); 
    sprintf('  ave. speed: %0.3f m/s', pathLen_cm/100/t_dur); 
    sprintf('  -----');    
    sprintf('  Output file: %s', vcFile_trial);
    sprintf('  Date analyzed: %s (%d days ago)', get_(S_dir, 'date'), nDaysAgo);
    sprintf('  version used: %s (%s)', handles.vcVer, handles.vcVer_date);
    sprintf('  -----');
    sprintf('  Video file: %s', get_(handles, 'vidFname'));   
    sprintf('    FPS: %0.3f', get_(handles, 'FPS'));       
    sprintf('  ADC file: %s', get_(handles, 'editADCfile', 'String')); 
    sprintf('  ADC_TS file: %s', get_(handles, 'editADCfileTs', 'String'));
};
% csMsg = [csMsg; get_(handles, 'csSettings')];
if nargout==0, disp(csMsg); end
end %func


%--------------------------------------------------------------------------
% 7/26/2018 JJJ: Copied from GUI.m
function [vcFile_full, S_dir] = fullpath_(vcFile)
[vcDir_, vcFile_, vcExt_] = fileparts(vcFile);
if isempty(vcDir_) 
    vcDir_ = pwd();
    vcFile_full = fullfile(vcDir_, vcFile);
else
    vcFile_full = vcFile;
end
if nargout>=2, S_dir = dir(vcFile_full); end
end %func


%--------------------------------------------------------------------------
function S_dir = file_dir_(vcFile_trial)
if exist_file_(vcFile_trial)
    S_dir = dir(vcFile_trial);    
else
    S_dir = [];
end
end %func

%--------------------------------------------------------------------------
function [mrXY_head, vrT] = get_traj_(handles)
P = load_settings_(handles);
Xs = filtPos(handles.XC, P.TRAJ_NFILT, 1);
Ys = filtPos(handles.YC, P.TRAJ_NFILT, 1);
mrXY_head = [Xs(:,2), Ys(:,2)];
vrT = get_(handles, 'TC');
end %func

        
%--------------------------------------------------------------------------
function commit_()
S_cfg = load_cfg_();
if exist_dir_('.git'), fprintf(2, 'Cannot commit from git repository\n'); return; end

% Delete previous files
S_warning = warning();
warning('off');
delete_empty_files_();
delete([S_cfg.vcDir_commit, '*']);
warning(S_warning);

% Copy files
copyfile_(S_cfg.csFiles_commit, S_cfg.vcDir_commit, '.');

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
vcVer = 'v0.2.3';
vcDate = '7/26/2018';
if nargout==0
    fprintf('%s (%s) installed\n', vcVer, vcDate);
    edit_('change_log.txt');
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
    '  vistrack download-sample';
    '    Download a sample video from Dropbox';
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
S_cfg = load_cfg_();
% delete_file_(get_(S_cfg, 'csFiles_delete'));

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
function hFig = figure_new_(vcTag, vcTitle)
if nargin<1, vcTag = ''; end
if nargin<2, vcTitle = ''; end

if ~isempty(vcTag)
    %remove prev tag duplication
    delete_multi_(findobj('Tag', vcTag, 'Type', 'Figure')); 
else
    hFig = figure('Tag', vcTag, 'Color', 'w', 'NumberTitle', 'off', 'Name', vcTitle);
end
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
for i=1:nargin()
    try
        delete(varargin{i});
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
% 7/19/18: Copied from jrc3.m
function val = get_set_(S, vcName, def_val)
% set a value if field does not exist (empty)

if isempty(S), S = get(0, 'UserData'); end
if isempty(S), val = def_val; return; end
if ~isstruct(S)
    val = []; 
    fprintf(2, 'get_set_: %s must be a struct\n', inputname(1));
    return;
end
val = get_(S, vcName);
if isempty(val), val = def_val; end
end %func


%--------------------------------------------------------------------------
% 7/19/18: Copied from jrc3.m
function varargout = get_(varargin)
% retrieve a field. if not exist then return empty
% [val1, val2] = get_(S, field1, field2, ...)
% val = get_(S, 'struct1', 'struct2', 'field');

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end
if nargout==1 && nargin > 2
    varargout{1} = get_recursive_(varargin{:}); return;
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
function out = get_recursive_(varargin)

% recursive get
out = [];
if nargin<2, return; end
S = varargin{1};
for iField = 2:nargin
    try
        out = S.(varargin{iField});
        if iField == nargin, return; end
        S = out;
    catch
        out = [];
    end
end % for
end %func


%--------------------------------------------------------------------------
% 7/19/2018
function P = load_settings_(handles)
% P = load_settings_()
% P = load_settings_(handles)

P = load_cfg_();
P.vcFile_settings = get_set_(P, 'vcFile_settings', 'settings_vistrack.m');
P.pixpercm = get_set_(P, 'pixpercm', 7.238);
P.angXaxis = get_set_(P, 'angXaxis', -0.946);

P_ = [];
try    
    csSettings = get(handles.editSettings, 'String');
    P_ = file2struct(csSettings);
catch
    P_ = file2struct(P.vcFile_settings);
end
P = struct_merge_(P, P_);
end %func


%--------------------------------------------------------------------------
% 7/19/2018 JJJ: Copied from jrc3.m
function P = struct_merge_(P, P1, csNames)
% Merge second struct to first one
% P = struct_merge_(P, P_append)
% P = struct_merge_(P, P_append, var_list) : only update list of variable names
if isempty(P), P=P1; return; end % P can be empty
if isempty(P1), return; end
if nargin<3, csNames = fieldnames(P1); end
if ischar(csNames), csNames = {csNames}; end

for iField = 1:numel(csNames)
    vcName_ = csNames{iField};
    if isfield(P1, vcName_), P.(vcName_) = P1.(vcName_); end
end
end %func


%--------------------------------------------------------------------------
function [mrPath, mrDur, S_trialset, cS_probe] = trialset_learningcurve_(vcFile_trialset)
% It loads the files
% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
% run S141106_LearningCurve_Control.m first cell

S_trialset = load_trialset_(vcFile_trialset);
[pixpercm, angXaxis] = struct_get_(S_trialset.P, 'pixpercm', 'angXaxis');
[tiImg, vcType_uniq, vcAnimal_uniq, viImg, csFiles_Track] = ...
    struct_get_(S_trialset, 'tiImg', 'vcType_uniq', 'vcAnimal_uniq', 'viImg', 'csFiles_Track');

[trDur, trPath, trFps] = deal(nan(size(tiImg)));
fprintf('Analyzing\n\t');
warning off;
t1 = tic;
cS_probe = {};
for iTrial = 1:numel(viImg)    
    try
        S_ = load(csFiles_Track{iTrial}, 'TC', 'XC', 'YC', 'xy0', 'vidFname', 'FPS', 'img0'); 
        iImg_ = viImg(iTrial);        
        if S_trialset.vlProbe(iTrial)
            cS_probe{end+1} = S_;
        else
            trPath(iImg_) = trial_pathlen_(S_, pixpercm, angXaxis);
            trDur(iImg_) = diff(S_.TC([1,end]));
        end
        trFps(iImg_) = get_set_(S_, 'FPS', nan);        
        fprintf('.');
    catch
        disp(csFiles_Track{iTrial});
    end
end %for
fprintf('\n\ttook %0.1fs\n', toc(t1));

% FPS integrity check
hFig = plot_trialset_img_(S_trialset, trFps); 
set(hFig, 'Name', sprintf('FPS: %s', vcFile_trialset));

% compact by removing nan. 
% date x session x animal (trPath,trDur) -> session x date x animal (trPath_,trDur_)
[nDates, nSessions, nAnimals] = size(tiImg);
[trPath_, trDur_] = deal(nan(nSessions, nDates, nAnimals));
for iAnimal = 1:size(tiImg,3)
    [mrPath1, mrDur1] = deal(trPath(:,:,iAnimal)', trDur(:,:,iAnimal)');
    vi1 = find(~isnan(mrPath1));
    vi2 = 1:numel(vi1);
    [mrPath2, mrDur2] = deal(nan(nSessions, nDates));
    [mrPath2(vi2), mrDur2(vi2)] = deal(mrPath1(vi1), mrDur1(vi1));    
    trPath_(:,:,iAnimal) = mrPath2;
    trDur_(:,:,iAnimal) = mrDur2;
end
[trPath_, trDur_] = deal(permute(trPath_,[1,3,2]), permute(trDur_,[1,3,2]));
[mrPath, mrDur] = deal(reshape(trPath_,[],nDates)/100, reshape(trDur_,[],nDates));
viCol = find(~any(isnan(mrPath)));
[mrPath, mrDur] = deal(mrPath(:,viCol), mrDur(:,viCol));

if nargout==0
    figure_new_('', ['Learning curve: ', vcFile_trialset]);
    subplot 211; errorbar_iqr_(mrPath); ylabel('Dist (m)'); grid on; xlabel('Session #');
    subplot 212; errorbar_iqr_(mrDur); ylabel('Duration (s)'); grid on; xlabel('Sesision #');
end
end %func


%--------------------------------------------------------------------------
function mr_ = errorbar_iqr_(mr)
mr_ = quantile_mr_(mr, [.25,.5,.75]);
errorbar(1:size(mr_,1), mr_(:,2), mr_(:,2)-mr_(:,1), mr_(:,3)-mr_(:,2));
set(gca, 'XLim', [.5, size(mr_,1)+.5]);
end %func


%--------------------------------------------------------------------------
function mr1 = quantile_mr_(mr, vrQ)
mr1 = zeros(numel(vrQ), size(mr,2), 'like', mr);
for i=1:size(mr,2)
    mr1(:,i) = quantile_(mr(:,i), vrQ);
end
mr1 = mr1';
end %func


%--------------------------------------------------------------------------
function [pathLen_cm, XH, YH, TC1] = trial_pathlen_(S_trial, pixpercm, angXaxis)
if nargin<2
    P = load_settings_();
    [pixpercm, angXaxis] = deal(P.pixpercm, P.angXaxis);
end %if

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


%--------------------------------------------------------------------------
function export_(handles)
assignWorkspace_(handles);
msgbox('"handles" struct assigned in workspace.');
disp(handles);
end %func


%--------------------------------------------------------------------------
% 7/20/2018 JJJ: list trialset files
function trialset_list_(vcFile_trialset)
S_trialset = load_trialset_(vcFile_trialset);
if isempty(S_trialset)
	errordlg(sprintf('%s does not exist', vcFile_trialset)); return; 
end
if ~exist_dir_(get_(S_trialset, 'vcDir'))
    errordlg(sprintf('vcDir=''%s''; does not exist', vcFile_trialset), vcFile_trialset); 
    return;
end
% S_trialset = load_trialset_(vcFile_trialset);
[tiImg, vcType_uniq, vcAnimal_uniq, csDir_trial, csFiles_Track] = ...
    struct_get_(S_trialset, 'tiImg', 'vcType_uniq', 'vcAnimal_uniq', 'csDir_trial', 'csFiles_Track');

% output
msgbox(S_trialset.csMsg, file_part_(vcFile_trialset));
disp_cs_(S_trialset.csMsg);
disp_cs_(S_trialset.csMsg2);

hFig = plot_trialset_img_(S_trialset, tiImg);
set(hFig, 'Name', sprintf('Integrity check: %s', vcFile_trialset));
end %func


%--------------------------------------------------------------------------
function hFig = plot_trialset_img_(S_trialset, tiImg, clim)
if nargin<2, tiImg = S_trialset.tiImg; end
if nargin<3, clim = [min(tiImg(:)), max(tiImg(:))]; end

hFig = figure_new_('', S_trialset.vcFile_trialset);
for iAnimal = 1:size(tiImg,3)
    subplot(1,size(tiImg,3),iAnimal);
    imagesc_(tiImg(:,:,iAnimal), clim);
    ylabel('Dates'); xlabel('Trials'); 
    title(sprintf('Animal %s', S_trialset.vcAnimal_uniq(iAnimal)));    
end %for
end %func


%--------------------------------------------------------------------------
function varargout = struct_get_(varargin)
% deal struct elements

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end

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
function S_trialset = load_trialset_(vcFile_trialset)
% return [] if vcFile_trialset does not exist or 
if ~exist_file_(vcFile_trialset), S_trialset = []; return; end
S_trialset = file2struct(vcFile_trialset);
P = load_settings_();
[csFiles_Track, csDir_trial] = find_files_(S_trialset.vcDir, '*_Track.mat');
if isempty(csFiles_Track), return; end
    
[csDataID, S_trialset_]  = get_dataid_(csFiles_Track);
S_trialset = struct_merge_(S_trialset, S_trialset_);
[vcAnimal_uniq, vnAnimal_uniq] = unique_(S_trialset.vcAnimal);
[viDate_uniq, vnDate_uniq] = unique_(S_trialset.viDate);
[vcType_uniq, vnType_uniq] = unique_(S_trialset.vcType);
[viTrial_uniq, vnTrial_uniq] = unique_(S_trialset.viTrial);
fh1_ = @(x,y,z)cell2mat(arrayfun(@(a,b)sprintf(z,a,b),x,y,'UniformOutput',0));
fh2_ = @(cs_)cell2mat(cellfun(@(vc_)sprintf('  %s\n',vc_),cs_,'UniformOutput',0));
csMsg = { ...
    sprintf('Trial types(#trials): %s\n', fh1_(vcType_uniq, vnType_uniq, '%c(%d), '));    
    sprintf('Animals(#trials): %s\n', fh1_(vcAnimal_uniq, vnAnimal_uniq, '%c(%d), '));
    sprintf('Dates(#trials):\n  %s\n', fh1_(viDate_uniq, vnDate_uniq, '%d(%d), '));
    sprintf('# Probe trials: %d', sum(S_trialset.vlProbe));
    sprintf('%s', fh2_(csFiles_Track(S_trialset.vlProbe)));
    sprintf('Figure color scheme: blue:no data, green:analyzed, yellow:probe trial');    
    sprintf('See the console output for further details');    
};

% image output
tiImg = zeros(max(viDate_uniq), max(viTrial_uniq), numel(vcAnimal_uniq));
viDate = S_trialset.viDate;
viAnimal = S_trialset.vcAnimal - min(S_trialset.vcAnimal) + 1;
viTrial = S_trialset.viTrial;
viImg = sub2ind(size(tiImg), viDate, viTrial, viAnimal);
tiImg(viImg) = 1;
tiImg(viImg(S_trialset.vlProbe)) = 2;

% find missing trials
[viDate_missing, viTrial_missing, viAnimal_missing] = ind2sub(size(tiImg), find(tiImg==0));
csDataID_missing = arrayfun(@(a,b,c)sprintf('%c%02d%c%d',vcType_uniq(1),a,b,c), ...
    viDate_missing, vcAnimal_uniq(viAnimal_missing)', viTrial_missing, ...
        'UniformOutput', 0);

fh3_ = @(cs)(cell2mat(cellfun(@(x)sprintf('  %s\n',x),cs, 'UniformOutput', 0)));
fh4_ = @(cs)(cell2mat(cellfun(@(x)sprintf('%s ',x),cs, 'UniformOutput', 0)));

% secondary message 
csMsg2 = { ...
    sprintf('\n[Folders]');
    fh3_(csDir_trial);
    sprintf('[Files]');
    fh3_(csFiles_Track);
    sprintf('[Probe trials]');
    fh2_(csFiles_Track(S_trialset.vlProbe));
    sprintf('[Missing trials (%d)]', numel(csDataID_missing));
    ['  ', fh4_(sort(csDataID_missing'))]
};

S_trialset = struct_add_(S_trialset, vcFile_trialset, P, ...
    csFiles_Track, csDir_trial, csMsg, csMsg2, ...
    tiImg, viDate, viTrial, viAnimal, viImg, ...
    vcAnimal_uniq, viDate_uniq, vcType_uniq, viTrial_uniq);
end %func


%--------------------------------------------------------------------------
function S = struct_add_(S, varargin)

for i=1:numel(varargin)
    S.(inputname(i+1)) = varargin{i};
end
end %func


%--------------------------------------------------------------------------
function imagesc_(mr, clim)
if nargin<2, clim = []; end
if isempty(clim)
    imagesc(mr, 'xdata', 1:size(mr,2), 'ydata', 1:size(mr,1));
else
    imagesc(mr, 'xdata', 1:size(mr,2), 'ydata', 1:size(mr,1), clim);
end
set(gca,'XTick', 1:size(mr,2));
set(gca,'YTick', 1:size(mr,1));
axis([.5, size(mr,2)+.5, .5, size(mr,1)+.5]);
grid on;
end %func


%--------------------------------------------------------------------------
function vc = file_part_(vc)
[~,a,b] = fileparts(vc);
vc = [a, b];
end %func


%--------------------------------------------------------------------------
function [vi_uniq, vn_uniq] = unique_(vi)
[vi_uniq, ~, vi_] = unique(vi);
vn_uniq = hist(vi_, 1:numel(vi_uniq));
end %func


%--------------------------------------------------------------------------
function [csDataID, S]  = get_dataid_(csFiles)
csDataID = cell(size(csFiles));
[viDate, viTrial] = deal(zeros(size(csFiles)));
vlProbe = false(size(csFiles));
[vcType, vcAnimal] = deal(repmat(' ', size(csFiles)));
for iFile=1:numel(csFiles)
    [~, vcFile_, ~] = fileparts(csFiles{iFile});
    vcDateID_ = strrep(vcFile_, '_Track', '');
    csDataID{iFile} = vcDateID_;    
    vcType(iFile) = vcDateID_(1);
    vcAnimal(iFile) = vcDateID_(4);
    viDate(iFile) = str2num(vcDateID_(2:3));
    viTrial(iFile) = str2num(vcDateID_(5));
    vlProbe(iFile) = numel(vcDateID_) > 5;
end %for
S = makeStruct_(vcType, viDate, vcAnimal, viTrial, vlProbe);
end %func


%--------------------------------------------------------------------------
% 7/20/18: Copied from jrc3.m
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name
S = struct();
for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
end %func


%--------------------------------------------------------------------------
function [csFiles, csDir] = find_files_(csDir, vcFile)
% consider using (dir('**/*.mat') for example instead of finddir

if ischar(csDir)
    if any(csDir=='*')
        csDir = find_dir_(csDir);
    else
        csDir = {csDir}; 
    end
end
csFiles = {};
for iDir=1:numel(csDir)
    vcDir_ = csDir{iDir};
    S_dir_ = dir(fullfile(vcDir_, vcFile));
    csFiles_ = cellfun(@(x)fullfile(vcDir_, x), {S_dir_.name}, 'UniformOutput', 0);
    csFiles = [csFiles, csFiles_];
end %for
end %func


%--------------------------------------------------------------------------
function csDir = find_dir_(vcDir)
% accepts if vcDir contains a wildcard
if ~any(vcDir=='*'), csDir = {vcDir}; return; end
[vcDir_, vcFile_, vcExt_] = fileparts(vcDir);
if ~isempty(vcExt_), csDir = {vcDir_}; return ;end

S_dir = dir(vcDir);
csDir = {S_dir.name};
csDir_ = csDir([S_dir.isdir]);
csDir = cellfun(@(x)fullfile(vcDir_, x), csDir_, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function vidobj = VideoReader_(vcFile_vid, nRetry)
if nargin<2, nRetry = []; end
if isempty(nRetry), nRetry = 3; end % number of frames can change
nThreads = 1; % disable parfor by setting it to 1. Parfor is slower

fprintf('Loading Video: %s\n', vcFile_vid); t1=tic;
cVidObj = cell(nRetry,1);
fParfor = is_parfor_(nThreads);
if fParfor
    try
        parfor iRetry = 1:nRetry
            [cVidObj{iRetry}, vnFrames(iRetry)] = load_vid_(vcFile_vid);
            fprintf('\t#%d: %d frames\n', iRetry, vnFrames(iRetry));
        end %for
    catch
        fParfor = 0;
    end
end
if ~fParfor
    for iRetry = 1:nRetry
        [cVidObj{iRetry}, vnFrames(iRetry)] = load_vid_(vcFile_vid);
        fprintf('\t#%d: %d frames\n', iRetry, vnFrames(iRetry));
    end %for
end
[NumberOfFrames, iMax] = max(vnFrames);
vidobj = cVidObj{iMax};

fprintf('\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [vidobj, nFrames] = load_vid_(vcFile_vid);
try
    vidobj = VideoReader(vcFile_vid);
    nFrames = vidobj.NumberOfFrames;
catch
    vidobj = [];
    nFrames = 0;    
end
end %func


%--------------------------------------------------------------------------
function fParfor = is_parfor_(nThreads)
if nargin<1, nThreads = []; end

if nThreads == 1
    fParfor = 0;
else
    fParfor = license('test', 'Distrib_Computing_Toolbox');
end
end %func


%--------------------------------------------------------------------------
% 11/5/17 JJJ: Created
function vc = dir_filesep_(vc)
% replace the file seperaation characters
if isempty(vc), return; end
vl = vc == '\' | vc == '/';
if any(vl), vc(vl) = filesep(); end
end %func


%--------------------------------------------------------------------------
function trialset_barplots_(vcFile_trialset)
% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
% run S141106_LearningCurve_Control.m first cell

[mrPath, mrDur, S_trialset] = trialset_learningcurve_(vcFile_trialset);

viEarly = get_(S_trialset, 'viEarly_trial');
viLate = get_(S_trialset, 'viLate_trial');
if isempty(viEarly) || isempty(viLate)
    msgbox('Set "viEarly_trial" and "viLate_trial" in .trialset file');
    return;
end

[vrPath_early, vrPath_late] = deal(mrPath(:,viEarly), mrPath(:,viLate));
[vrDur_early, vrDur_late] = deal(mrDur(:,viEarly), mrDur(:,viLate));
[vrSpeed_early, vrSpeed_late] = deal(vrPath_early./vrDur_early, vrPath_late./vrDur_late);
quantLim = get_set_(S_trialset, 'quantLim', [1/8, 7/8]);
[vrPath_early, vrPath_late, vrDur_early, vrDur_late, vrSpeed_early, vrSpeed_late] = ...
    trim_quantile_(vrPath_early, vrPath_late, vrDur_early, vrDur_late, vrSpeed_early, vrSpeed_late, quantLim);

figure_new_('', ['Early vs Late: ', vcFile_trialset]);
subplot 131;
bar_mean_sd_({vrPath_early, vrPath_late}, {'Early', 'Late'}, 'Pathlen (m)');
subplot 132;
bar_mean_sd_({vrDur_early, vrDur_late}, {'Early', 'Late'}, 'Duration (s)');
subplot 133;
bar_mean_sd_({vrSpeed_early, vrSpeed_late}, {'Early', 'Late'}, 'Speed (m/s)');

msgbox(sprintf('Early Sessions: %s\nLate Sessions: %s', sprintf('%d ', viEarly), sprintf('%d ', viLate)));
end %func


%--------------------------------------------------------------------------
function varargout = bar_mean_sd_(cvr, csXLabel, vcYLabel)
if nargin<2, csXLabel = {}; end
if nargin<3, vcYLabel = ''; end
if isempty(csXLabel), csXLabel = 1:numel(cvr); end

vrMean = cellfun(@(x)nanmean(x(:)), cvr);
vrSd = cellfun(@(x)nanstd(x(:)), cvr);
vrX = 1:numel(cvr);

errorbar(vrX, vrMean, [], vrSd, 'k', 'LineStyle', 'none'); 
hold on; grid on;
h = bar(vrX, vrMean);
set(h, 'EdgeColor', 'None');
set(gca, 'XTick', vrX, 'XTickLabel', csXLabel, 'XLim', vrX([1,end]) + [-.5, .5]);
ylabel(vcYLabel);

[h,pa]=ttest2(cvr{1},cvr{2});
fprintf('%s: E vs L, p=%f\n', vcYLabel, pa);
end %func


%--------------------------------------------------------------------------
function varargout = trim_quantile_(varargin)
qlim = varargin{end};
for iArg = 1:nargout
    vr_ = varargin{iArg};
    varargout{iArg} = quantFilt_(vr_(:), qlim);
end %for
end %func


%--------------------------------------------------------------------------
function vr = quantFilt_(vr, quantLim)
qlim = quantile(vr(:), quantLim);
vr = vr(vr >= qlim(1) & vr < qlim(end));
end %func


%--------------------------------------------------------------------------
% Display list of toolbox and files needed
% 7/26/17 JJJ: Code cleanup and test
function [fList, pList] = disp_dependencies_(vcFile)
if nargin<1, vcFile = []; end
if isempty(vcFile), vcFile = mfilename(); end

[fList,pList] = matlab.codetools.requiredFilesAndProducts(vcFile);
if nargout==0
    disp('Required toolbox:');
    disp({pList.Name}');
    disp('Required files:');
    disp(fList');
end
end % func


%--------------------------------------------------------------------------
function download_sample_()
S_cfg = load_cfg_();
csLink = get_(S_cfg, 'csLink_sample');
if isempty(csLink), fprintf(2, 'Sample video does not exist\n'); return; end

t1 = tic;
fprintf('Downloading sample files. This can take up to several minutes.\n');
vlSuccess = download_files_(csLink);
fprintf('\t%d/%d files downloaded. Took %0.1fs\n', ...
    sum(vlSuccess), numel(vlSuccess), toc(t1));
end %func


%--------------------------------------------------------------------------
function vlSuccess = download_files_(csLink, csDest)
% download file from the web
nRetry = 5;

if nargin<2, csDest = link2file_(csLink); end
vlSuccess = false(size(csLink));
for iFile=1:numel(csLink)    
    for iRetry = 1:nRetry
        try
            % download from list of files    
            fprintf('\tDownloading %s: ', csLink{iFile});
            vcFile_out1 = websave(csDest{iFile}, csLink{iFile});
            fprintf('saved to %s\n', vcFile_out1);
            vlSuccess(iFile) = 1;
            break;
        catch
            fprintf('\tRetrying %d/%d\n', iRetry, nRetry);
            if iRetry==nRetry
                fprintf(2, '\n\tDownload failed. Please download manually from the link below.\n');
                fprintf(2, '\t%s\n', csLink{iFile});
            end
        end
    end
end %for
end %func


%--------------------------------------------------------------------------
function csFile = link2file_(csLink)
csFile = cell(size(csLink));
for i=1:numel(csLink)        
    vcFile1 = csLink{i};
    iBegin = find(vcFile1=='/', 1, 'last'); % strip ?    
    if ~isempty(iBegin), vcFile1 = vcFile1(iBegin+1:end); end

    iEnd = find(vcFile1=='?', 1, 'last'); % strip ?
    if ~isempty(iEnd), vcFile1 = vcFile1(1:iEnd-1); end
    csFile{i} = vcFile1;
end
end %func


%--------------------------------------------------------------------------
% 7/25/2018 JJJ: Wait for file to get deleted
function delete_file_(csFiles)
if isempty(csFiles), return; end
if ischar(csFiles), csFiles = {csFiles}; end
nRetry = 5;

for iRetry = 1:nRetry
    for iFile = 1:numel(csFiles)
        if ~exist_file_(csFiles{iFile}), continue; end
        delete_(csFiles{iFile});
    end
end
for i=1:nRetry, pause(.2); end % wait for file deletion
end %func


%--------------------------------------------------------------------------
% 7/25/2018 JJJ: Wait for file to get deleted
function S_cfg = load_cfg_()
try
    S_cfg = file2struct('default.cfg');    
catch
    S_cfg = struct(); % return an empty struct
end

% default field
S_cfg.vcDir_commit = get_set_(S_cfg, 'vcDir_commit', 'D:\Dropbox\Git\vistrack\'); 
S_cfg.csFiles_commit = get_set_(S_cfg, 'csFiles_commit', {'*.m', 'GUI.fig', 'change_log.txt', 'readme.txt', 'example.trialset', 'default.cfg'});
S_cfg.csFiles_delete = get_set_(S_cfg, 'csFiles_delete', {'settings_vistrack.m', 'example.trialset', 'R12A2_Track.mat'});
S_cfg.quantLim = get_set_(S_cfg, 'quantLim', [1/8, 7/8]);
end %func

