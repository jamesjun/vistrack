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
    case 'load-cfg', varargout{1} = load_cfg_();
    case 'clear-cache', clear_cache_();
    case 'loadvid-preview', varargout{1} = loadvid_preview_(vcArg1, vcArg2);
        
    case 'trial-sync', varargout{1} = trial_sync_(vcArg1);
    case 'cam2adc-sync', varargout{1} = cam2adc_sync_(vcArg1, vcArg2);
    case 'adc2cam-sync', varargout{1} = adc2cam_sync_(vcArg1, vcArg2);
    case 'trial-visitcount', trial_timemap_(vcArg1);
    case 'trial-fixsync', varargout{1} = trial_fixsync_(vcArg1);
    case 'trial-save', varargout{1} = trial_save_(vcArg1);
        
    case 'trialset-list', trialset_list_(vcArg1);    
    case 'trialset-learningcurve', trialset_learningcurve_(vcArg1);
    case 'trialset-barplots', trialset_barplots_(vcArg1);
%     case 'trialset-probe', trialset_probe_(vcArg1);
    case 'trialset-exportcsv', trialset_exportcsv_(vcArg1);
    case 'trialset-checkfps', trialset_checkfps_(vcArg1);
    case 'trialset-coordinates', trialset_coordinates_(vcArg1);
    case 'trialset-fixfps', trialset_fixfps_(vcArg1);
    case 'trialset-import-track', trialset_import_track_(vcArg1);
    case {'trialset-googlesheet', 'googlesheet'}, trialset_googlesheet_(vcArg1);

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
vcVer = 'v0.3.7';
vcDate = '8/29/2018';
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
        system(['start ', vcPage]);
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
function flag = exist_file_(vcFile, fVerbose)
if nargin<2, fVerbose = 0; end
if ~ischar(vcFile), flag = 0; return; end
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
fOverwrite = 1;

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
        if fOverwrite
            code = system('git fetch --all'); 
            code = system('git reset --hard origin/master'); 
        else
            code = system('git pull'); % do not overwrite existing changes
        end
    else
        code = system('git fetch --all');         
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
% 8/7/2018 JJJ
function flag = exist_dir_(vcDir)
if isempty(vcDir)
    flag = 0;
else
    S_dir = dir(vcDir);
    if isempty(S_dir)
        flag = 0;
    else
        flag = sum([S_dir.isdir]) > 0;
    end
%     flag = exist(vcDir, 'dir') == 7;
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
function close_(varargin)
for i=1:nargin
    v_ = varargin{i};
    try 
        if iscell(v_)
            close_(v_{:});
        elseif numel(v_)>1
            close_(v_);
        elseif numel(v_) == 1
            close(v_); 
        end
    catch 
        ;
    end
end
end %func


%--------------------------------------------------------------------------
function hFig = figure_new_(vcTag, vcTitle, vrPos)
if nargin<1, vcTag = ''; end
if nargin<2, vcTitle = ''; end
if nargin<3, vrPos = []; end

if ~isempty(vcTag)
    %remove prev tag duplication
    delete_multi_(findobj('Tag', vcTag, 'Type', 'Figure')); 
end
hFig = figure('Tag', vcTag, 'Color', 'w', 'NumberTitle', 'off', 'Name', vcTitle);
if ~isempty(vrPos), resize_figure_(hFig, vrPos); drawnow; end
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
        v_ = varargin{i};
        if iscell(v_)
            delete_(v_{:});
        elseif numel(v_) > 1
            for i=1:numel(v_), delete_(v_(i)); end
        elseif numel(v_) == 1
            delete(v_);
        end
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
if nargin<1, handles = []; end
P = load_cfg_();
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
function [mrPath, mrDur, S_trialset, cS_trial] = trialset_learningcurve_(vcFile_trialset)
% It loads the files
% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
% run S141106_LearningCurve_Control.m first cell

S_trialset = load_trialset_(vcFile_trialset);
[pixpercm, angXaxis] = struct_get_(S_trialset.P, 'pixpercm', 'angXaxis');
[tiImg, vcType_uniq, vcAnimal_uniq, viImg, csFiles_Track] = ...
    struct_get_(S_trialset, 'tiImg', 'vcType_uniq', 'vcAnimal_uniq', 'viImg', 'csFiles_Track');

hMsg = msgbox('Analyzing... (This closes automatically)');
[trDur, trPath, trFps] = deal(nan(size(tiImg)));
fprintf('Analyzing\n\t');
warning off;
t1 = tic;
cS_trial = {};
csField_load = setdiff(S_trialset.P.csFields, {'MOV', 'ADC','img1','img00'}); % do not load 'MOV' field since it's big and unused
for iTrial = 1:numel(viImg)    
    try
        S_ = load(csFiles_Track{iTrial}, csField_load{:}); 
        S_.vcFile_Track = csFiles_Track{iTrial};
        iImg_ = viImg(iTrial);                
        cS_trial{end+1} = S_;
        if ~S_trialset.vlProbe(iTrial)
            trPath(iImg_) = trial_pathlen_(S_, pixpercm, angXaxis);
            trDur(iImg_) = diff(S_.TC([1,end]));
        end
        trFps(iImg_) = get_set_(S_, 'FPS', nan);        
        fprintf('.');
    catch
        fprintf(2, '\n\tExport error: %s\n\t%s\n', csFiles_Track{iTrial}, lasterr());
    end
end %for
fprintf('\n\ttook %0.1fs\n', toc(t1));
close_(hMsg);

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
[trPath_, trDur_] = deal(permute(trPath_,[1,3,2]), permute(trDur_,[1,3,2])); % nSessions x nAnimals x nDate 
[mrPath, mrDur] = deal(reshape(trPath_,[],nDates)/100, reshape(trDur_,[],nDates));
viCol = find(~any(isnan(mrPath)));
[mrPath, mrDur] = deal(mrPath(:,viCol), mrDur(:,viCol));

% Export to csv
vcAnimals = cell2mat(S_trialset.csAnimals);
vcFile_path = subsFileExt_(vcFile_trialset, sprintf('_pathlen_%s.csv', vcAnimals));
vcFile_dur = subsFileExt_(vcFile_trialset, sprintf('_duration_%s.csv', vcAnimals));
csvwrite_(vcFile_path, mrPath', 'Learning-curve path-length (m), [sessions, trials]');
csvwrite_(vcFile_dur, mrDur', 'Learning-curve: duration (s), [sessions, trials]');

if nargout==0
    % FPS integrity check
    hFig = plot_trialset_img_(S_trialset, trFps); 
    set(hFig, 'Name', sprintf('FPS: %s', vcFile_trialset));
    
    % Plot learning curve
    figure_new_('', ['Learning curve: ', vcFile_trialset, '; Animals:', cell2mat(S_trialset.csAnimals)]);
    subplot 211; errorbar_iqr_(mrPath); ylabel('Dist (m)'); grid on; xlabel('Session #');
    subplot 212; errorbar_iqr_(mrDur); ylabel('Duration (s)'); grid on; xlabel('Sesision #');    
end
end %func


%--------------------------------------------------------------------------
function vcMsg = csvwrite_(vcFile, var, vcMsg)
if nargin<3, vcVar = inputname(2); end
try
    csvwrite(vcFile, var);
    vcMsg = sprintf('"%s" wrote to %s', vcMsg, vcFile);
catch
    vcMsg = sprintf('Failed to write "s" to %s', vcMsg, vcFile);
end
if nargout==0, disp(vcMsg); end
end %func


%--------------------------------------------------------------------------
function [S_trialset, trFps] = trialset_checkfps_(vcFile_trialset)
% It loads the files
% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
% run S141106_LearningCurve_Control.m first cell
% fFix_sync = 0;

S_trialset = load_trialset_(vcFile_trialset);
% [pixpercm, angXaxis] = struct_get_(S_trialset.P, 'pixpercm', 'angXaxis');
[tiImg, vcType_uniq, vcAnimal_uniq, viImg, csFiles_Track] = ...
    struct_get_(S_trialset, 'tiImg', 'vcType_uniq', 'vcAnimal_uniq', 'viImg', 'csFiles_Track');
warning off;
hMsg = msgbox('Analyzing... (This closes automatically)');
t1=tic;
trFps = nan(size(tiImg));
for iTrial = 1:numel(viImg)    
    try
        S_ = load(csFiles_Track{iTrial}, 'TC', 'XC', 'YC', 'xy0', 'vidFname', 'FPS', 'img0'); 
        if isempty(get_(S_, 'FPS')) || isempty(get_(S_, 'TC')), error('FPS or TC not found'); end
        S_.vcFile_Track = csFiles_Track{iTrial};
%         if fFix_sync, S_ = trial_fixsync_(S_, 0); end
        iImg_ = viImg(iTrial);        
        trFps(iImg_) = get_set_(S_, 'FPS', nan);                
        fprintf('.');
    catch
        fprintf('\n\t%s: %s\n', csFiles_Track{iTrial}, lasterr());
    end
end %for
fprintf('\n\ttook %0.1fs\n', toc(t1));
close_(hMsg);

if nargout==0
    hFig = plot_trialset_img_(S_trialset, trFps); 
    set(hFig, 'Name', sprintf('FPS: %s', vcFile_trialset));
end
msgbox('Click to display file locations (copies to the clipboard)');
% open a FPS fix tool 

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
    mr1(:,i) = quantile(mr(:,i), vrQ);
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
[dataID, fishID] = trial_id_(S_trial);
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
% export heading angles to CVS file
[vcFile_cvs, mrTraj, vcMsg_cvs, csFormat] = trial2csv_(handles, [], 0);
assignWorkspace_(mrTraj);
P = load_cfg_();
csMsg = {'"handles" struct and "mrTraj" assigned to the Workspace.', vcMsg_cvs};
csMsg = [csMsg, csFormat(:)'];
msgbox_(csMsg);
end %func


%--------------------------------------------------------------------------
% 10/15/2018 JJJ: cvs renamed to csv
function [vcFile_cvs, mrTraj, vcMsg, csFormat] = trial2csv_(S_trial, P, fPlot)
if nargin<2, P = []; end
if nargin<3, fPlot = 0; end
if isempty(P), P = load_settings_(); end

[vcDir_, ~, ~] = fileparts(S_trial.vidFname);
if exist_file_(get_(S_trial, 'vcFile_Track'))
    vcFile_cvs = subsFileExt_(S_trial.vcFile_Track, '.csv');
elseif exist_dir_(vcDir_)
    vcFile_cvs = subsFileExt_(S_trial.vidFname, '_Track.csv');
else
    vcFile_Track = get_(S_trial.editResultFile, 'String');
    vcFile_cvs = strrep(vcFile_Track, '.mat', '.csv');
end

P1 = setfield(P, 'xy0', S_trial.xy0);
mrTraj = resample_trial_(S_trial, P1);
csvwrite(vcFile_cvs, mrTraj);
vcMsg = sprintf('Trajectory exported to %s\n', vcFile_cvs);

% Export shape
if isfield(S_trial, 'mrPos_shape')
    vcFile_shapes = strrep(vcFile_cvs, '_Track.csv', '_shapes.csv');
    cm_per_grid = get_set_(P, 'cm_per_grid', 5);
    mrPos_shape_meter = S_trial.mrPos_shape;
    mrPos_shape_meter(:,1:2) = mrPos_shape_meter(:,1:2) * cm_per_grid / 100;
    csvwrite(vcFile_shapes, mrPos_shape_meter);
    vcMsg = sprintf('%sShapes exported to %s\n', vcMsg, vcFile_shapes);

    vcFile_relations = strrep(vcFile_cvs, '_Track.csv', '_relations.csv');    
    mrRelations = calc_relations_(mrTraj, mrPos_shape_meter, P1);
    csvwrite(vcFile_relations, mrRelations);
    vcMsg = sprintf('%srelations exported to %s\n', vcMsg, vcFile_relations);
else
    fprintf(2, '%s: Shape positions field does not exist.\n', vcFile_cvs);
end

csShapes = get_(P, 'csShapes');
csFormat = {...
    '_Track.csv files:', 
    '  Columns: T(s), X(m), Y(m), A(deg), R(Hz), D(m), V(m/s), S(m):',
    '    T: camera frame time',
    '    X: x coordinate of the head tip @ grid frame of reference',
    '    Y: y coordinate of the head tip @ grid frame of reference',
    '    A: head orientation',
    '    R: EOD rate',
    '    D: Distance per EOD pulse (=1/sampling_density)',
    '    V: Head speed (m/s, signed)',
    '    S: Distance per Escan (=1/escan_density)',    
    '_shapes.csv files:', 
    '  Columns: X(m), Y(m), A(deg):',
    '    X(m): x coordinate of the shape center @ grid frame of reference',
    '    Y(m): y coordinate of the shape center @ grid frame of reference',
    '    A(deg): Shape orientation',
    sprintf('  Rows: %s', sprintf('"%s", ', csShapes{:})),     
    '_relations.csv files:',
    sprintf('  Columns: T(s), D_F(m), A_E(deg), %s', sprintf('L_"%s"(bool), ', csShapes{:})),
    '    T: camera frame time',
    '    D_F: distance to the food', 
    '    A_E: heading angle error (food_vec - head_vec, 0..90 deg)',    
    '    L_"x": Is shape "x" adjacent to the head position? 0:no, 1:yes'};    

if fPlot
    hFig = figure_new_('', vcFile_cvs);
    imshow(S_trial.img0); hold on;
    resize_figure_(hFig, [0,0,.5,1]);
    plot_chevron_(S_trial.XC(:,2:3), S_trial.YC(:,2:3));
end
if nargout==0
    fprintf('%s', vcMsg); 
    disp_cs_(csFormat);
end
end %func


%--------------------------------------------------------------------------
function [mrTXYARDVS_rs, P1] = resample_trial_(S_trial, P)
% Output
% -----
% mrTXYARDV_rs: 
%   Time(s), X-pos(m), Y-pos(m), Angle(deg), Sampling Rate(Hz),
%   Dist/pulse(m), Velocity(m/s), Dist/E-Scan (m)

% smooth the trajectory
% if fFilter
fh_filt = @(x)filtPos(x, P.TRAJ_NFILT, 1);
P1 = setfield(P, 'xy0', S_trial.xy0);
sRateHz_rs = get_set_(P, 'sRateHz_resample', 100);
vrT_rs = (S_trial.TC(1):1/sRateHz_rs:S_trial.TC(end))';
fh_interp = @(x)interp1(S_trial.TC(:), x, vrT_rs);

mrXY_pix_2 = [fh_filt(S_trial.XC(:,2)), fh_filt(S_trial.YC(:,2))];
mrXY_pix_3 = [fh_filt(S_trial.XC(:,3)), fh_filt(S_trial.YC(:,3))];
mrXY_m_2_rs = fh_interp(pix2cm_(mrXY_pix_2, P1) / 100);
mrXY_m_3_rs = fh_interp(pix2cm_(mrXY_pix_3, P1) / 100);
vrA_rs = fh_interp(pix2cm_deg_(S_trial.AC(:,2), P1));

% add EOD and sampling density
[vtEodr, vrEodr] = getEodr_adc_(S_trial, P);
vrR_rs = interp1(vtEodr, vrEodr, vrT_rs);

% Compute velocity
mrXY_23_rs = mrXY_m_2_rs - mrXY_m_3_rs;
[VX, VY] = deal(diff3_(mrXY_m_2_rs(:,1)), diff3_(mrXY_m_2_rs(:,2)));
vrV_rs = hypot(VX, VY) .* sign(mrXY_23_rs(:,1).*VX + mrXY_23_rs(:,2).*VY) * sRateHz_rs;

% Count sampling density
vrLr = cumsum(hypot(VX, VY));
vtEodr_ = vtEodr(vtEodr>=vrT_rs(1) & vtEodr <= vrT_rs(end));
vrD = diff3_(interp1(vrT_rs, vrLr, vtEodr_, 'spline')); % distance between EOD
vrD_rs = interp1(vtEodr_, vrD, vrT_rs);
vrD_rs(isnan(vrD_rs)) = 0;

% calc ESCAN rate
viDs = findDIsac(diff3_(diff3_(vtEodr))); 
vtEscan = vtEodr(viDs);
vrDs = diff3_(interp1(vrT_rs, vrLr, vtEscan, 'spline')); % distance between EOD
vrS_rs = interp1(vtEscan, vrDs, vrT_rs, 'spline');
vrS_rs(isnan(vrS_rs)) = 0;

% get EOD timestamps
mrTXYARDVS_rs = [vrT_rs, mrXY_m_2_rs, vrA_rs(:), vrR_rs(:), vrD_rs(:), vrV_rs(:), vrS_rs(:)];
% figure; quiver(mrXY_m_2_rs(:,1),mrXY_m_2_rs(:,2), VX, VY, 2, 'r.')
end %func


%--------------------------------------------------------------------------
function data = diff3_(data, dt) % three point diff
%data = differentiate5(data, dt)
%data: timeseries to differentiate
%dt: time step (default of 1)
% http://en.wikipedia.org/wiki/Numerical_differentiation
dimm = size(data);
data=data(:)';
data = filter([1 0 -1], 2, data);
data = data(3:end);
data = [data(1), data, data(end)];

if nargin > 1
    data = data / dt;
end

if dimm(1) == 1 %row vector
    data=data(:)';
else
    data=data(:); %col vector
end
end %func


%--------------------------------------------------------------------------
function mrRelations = calc_relations_(mrTraj, mrPos_shape_meter, P)
if nargin<3, P = []; end
if isempty(P), P = load_settings_(); end

[T, mrXY_h, A_H] = deal(mrTraj(:,1), mrTraj(:,2:3), mrTraj(:,4));
vrXY_food = mrPos_shape_meter(end,1:2);
mrV_F = bsxfun(@minus, vrXY_food, mrXY_h);
[A_F, D_F] = cart2pol_(mrV_F(:,1), mrV_F(:,2));
% A_E = min(mod(A_F-A_H, 180), mod(A_H-A_F, 180));
A_E = mod(A_F-A_H+90, 180) - 90;

% determine shape mask
dist_cut = get_set_(P, 'dist_cm_shapes', 3) / 100; % in meters
nShapes = size(mrPos_shape_meter,1);
mlL_shapes = false(numel(T), nShapes);
for iShape = 1:nShapes
    vcShape = strtok(P.csShapes{iShape}, ' ');
    xya_ = mrPos_shape_meter(iShape,:);
    len_ = P.vrShapes(iShape)/100;
    [mrXY_poly_, fCircle] = get_polygon_(vcShape, xya_(1:2), len_, xya_(3));
    if fCircle
        vrD_ = hypot(xya_(1)-mrXY_h(:,1), xya_(2)-mrXY_h(:,2)) - len_/2; 
    else %polygon
        vrD_ = nearest_perimeter_(mrXY_poly_/100, mrXY_h); % convert to meter
    end
    mlL_shapes(:,iShape) = vrD_ <= dist_cut;
end

mrRelations = [T, D_F, A_E, mlL_shapes];
end %func


%--------------------------------------------------------------------------
function vrD_ = nearest_perimeter_(mrXY_p, mrXY_h)
% mrXY_p: polygon vertices
nInterp = 100;

mrXY_ = [mrXY_p; mrXY_p(1,:)]; % wrap
mrXY_int = interp1(1:size(mrXY_,1), mrXY_, 1:1/nInterp:size(mrXY_,1));
vrD_ = min(pdist2(mrXY_int, mrXY_h))';

if nargout==0
    figure; hold on;
    plot(mrXY_int(:,1), mrXY_int(:,2), 'b.-');
    plot(mrXY_h(:,1), mrXY_h(:,2), 'r.-');
end
end %func


%--------------------------------------------------------------------------
function [th_deg, r] = cart2pol_(x,y)
% th: degrees
th_deg = atan2(y,x).*(180/pi);
r = hypot(x,y);
end %func


%--------------------------------------------------------------------------
function mrA1 = pix2cm_deg_(mrA, P1)
angXaxis = get_set_(P1, 'angXaxis', 0);
mrA1 = mod(-mrA - angXaxis,360);
end


%--------------------------------------------------------------------------
function plot_chevron_(mrX, mrY)
nFrames = size(mrX,1);
hold on;
for iFrame=1:nFrames
    plotChevron(mrX(iFrame,:), mrY(iFrame,:), [], 90, .3);
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: added '.' if dir is empty
% 7/31/17 JJJ: Substitute file extension
function varargout = subsFileExt_(vcFile, varargin)
% Substitute the extension part of the file
% [out1, out2, ..] = subsFileExt_(filename, ext1, ext2, ...)

[vcDir_, vcFile_, ~] = fileparts(vcFile);
if isempty(vcDir_), vcDir_ = '.'; end
for i=1:numel(varargin)
    vcExt_ = varargin{i};    
    varargout{i} = [vcDir_, filesep(), vcFile_, vcExt_];
end
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
csFiles_Track = get_(S_trialset, 'csFiles_Track');
if isempty(csFiles_Track)
    errordlg(sprintf('No _Track.mat files are found in "%s".', vcFile_trialset), vcFile_trialset); 
    return;
end
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
function [hFig, vhImg] = plot_trialset_img_(S_trialset, tiImg, clim)
% make a click callback and show video location

if nargin<2, tiImg = S_trialset.tiImg; end
if nargin<3, clim = [min(tiImg(:)), max(tiImg(:))]; end

vhImg = zeros(size(tiImg,3), 1);
hFig = figure_new_('FigOverview', S_trialset.vcFile_trialset, [.5,0,.5,1]);
set0_(S_trialset);
% set0_(S_trialset);
% vhAx = zeros(size(tiImg,3), 1);
for iAnimal = 1:size(tiImg,3)
    hAx_ = subplot(1,size(tiImg,3),iAnimal);
    hImg_ = imagesc_(tiImg(:,:,iAnimal), clim);    
    vhImg(iAnimal) = hImg_;
    ylabel('Dates'); xlabel('Trials'); 
    title(sprintf('Animal %s', S_trialset.vcAnimal_uniq(iAnimal)));
%     hAx_.UserData = [];
    hImg_.ButtonDownFcn = @(h,e)button_FigOverview_(h,e,iAnimal);
end %for
end %func


%--------------------------------------------------------------------------
function button_FigOverview_(hImg, event, iAnimal)
% S_axes = get(hAxes, 'UserData');
xy = get(hImg.Parent, 'CurrentPoint');
xy = round(xy(1,1:2));

[iSession, iTrial, cAnimal] = deal(xy(2), xy(1), 'A' + iAnimal - 1);

fprintf('Session:%d, Trial:%d, Animal:%c\n', iSession, iTrial, cAnimal);
S_trialset = get0_('S_trialset');
vcVidExt = get_set_(S_trialset.P, 'vcVidExt', '.wmv');
% vcFormat = sprintf('*%02d%c%d.%s$', iSession, cAnimal, iTrial, vcVidExt)
vcFormat = sprintf('%02d%c%d(\\w*)_Track.mat', iSession, cAnimal, iTrial);
cs = cellfun(@(x)regexpi(x, vcFormat, 'match'), S_trialset.csFiles_Track, 'UniformOutput', 0);
iFind = find(~cellfun(@isempty, cs));
if ~isempty(iFind)
    vcFile_Track = S_trialset.csFiles_Track{iFind};
    vcFile_vid = strrep(vcFile_Track, '_Track.mat', vcVidExt);    
    switch event.Button
        case 1
            clipboard('copy', vcFile_Track);
            fprintf('\t%s (copied)\n\t%s\n', vcFile_Track, vcFile_vid);            
        case 3
            clipboard('copy', vcFile_vid);
            fprintf('\t%s\n\t%s (copied)\n', vcFile_Track, vcFile_vid);
    end
end
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
P = load_settings_();
S_trialset = file2struct(vcFile_trialset);
[csFiles_Track, csDir_trial] = find_files_(S_trialset.vcDir, '*_Track.mat');
if isempty(csFiles_Track), S_trialset.P=P; return; end
    
[csDataID, S_trialset_, csFiles_Track]  = get_dataid_(csFiles_Track, get_(S_trialset, 'csAnimals'));
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
    viDate_missing, toVec_(vcAnimal_uniq(viAnimal_missing)), viTrial_missing, ...
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

S_trialset = struct_add_(S_trialset, P, vcFile_trialset, ...
    csFiles_Track, csDir_trial, csMsg, csMsg2, ...
    tiImg, viDate, viTrial, viAnimal, viImg, ...
    vcAnimal_uniq, viDate_uniq, vcType_uniq, viTrial_uniq);
end %func


%--------------------------------------------------------------------------
function vr = toVec_(vr)
vr = vr(:);
end %func


%--------------------------------------------------------------------------
function vr = toRow_(vr)
vr = vr(:)';
end %func


%--------------------------------------------------------------------------
function S = struct_add_(S, varargin)

for i=1:numel(varargin)
    S.(inputname(i+1)) = varargin{i};
end
end %func


%--------------------------------------------------------------------------
function hImg = imagesc_(mr, clim)
if nargin<2, clim = []; end
if isempty(clim)
    hImg = imagesc(mr, 'xdata', 1:size(mr,2), 'ydata', 1:size(mr,1));
else
    hImg = imagesc(mr, 'xdata', 1:size(mr,2), 'ydata', 1:size(mr,1), clim);
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
function [csDataID, S, csFiles]  = get_dataid_(csFiles, csAnimals)
if nargin<2, csAnimals = {}; end


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

% Filter by animals
if ~isempty(csAnimals)
    vcAnimal_plot = cell2mat(csAnimals);
    viKeep = find(ismember(vcAnimal, vcAnimal_plot));
    [vcType, viDate, vcAnimal, viTrial, vlProbe, csFiles] = ...
        deal(vcType(viKeep), viDate(viKeep), vcAnimal(viKeep), viTrial(viKeep), vlProbe(viKeep), csFiles(viKeep));
end

S = makeStruct_(vcType, viDate, vcAnimal, viTrial, vlProbe, csFiles);
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

fprintf('Opening Video: %s\n', vcFile_vid); t1=tic;
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

% [mrPath, mrDur, S_trialset, cS_trial] = trialset_learningcurve_(vcFile_trialset);
[cS_trial, S_trialset, mrPath, mrDur] = loadShapes_trialset_(vcFile_trialset);

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
vcAnimal_use = cell2mat(S_trialset.csAnimals);

figure_new_('', ['Early vs Late: ', vcFile_trialset, '; Animals: ', vcAnimal_use]);
subplot 131;
bar_mean_sd_({vrPath_early, vrPath_late}, {'Early', 'Late'}, 'Pathlen (m)');
subplot 132;
bar_mean_sd_({vrDur_early, vrDur_late}, {'Early', 'Late'}, 'Duration (s)');
subplot 133;
bar_mean_sd_({vrSpeed_early, vrSpeed_late}, {'Early', 'Late'}, 'Speed (m/s)');
msgbox(sprintf('Early Sessions: %s\nLate Sessions: %s', sprintf('%d ', viEarly), sprintf('%d ', viLate)));

% Plot probe trials
S_shape = pool_probe_trialset_(S_trialset, cS_trial);

vcFigName = sprintf('%s; Animals: %s; Probe trials', S_trialset.vcFile_trialset, cell2mat(S_trialset.csAnimals));
hFig = figure_new_('', vcFigName, [0 .5 .5 .5]); 
viShapes = 1:6;
subplot 241; bar(1./S_shape.mrDRVS_shape(1,viShapes)); ylabel('Sampling density'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
subplot 242; bar(S_shape.mrDRVS_shape(2,viShapes)); ylabel('Sampling Rate (Hz)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
subplot 243; bar(S_shape.mrDRVS_shape(3,viShapes)); ylabel('Speed (m/s)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
subplot 244; bar(1./S_shape.mrDRVS_shape(4,viShapes)); ylabel('EScan density'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;        
subplot 245; bar(S_shape.vnVisit_shape(viShapes)); ylabel('# visits'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
subplot 246; bar(S_shape.vtVisit_shape(viShapes)); ylabel('t visit (s)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
subplot 247; bar(S_shape.vtVisit_shape(viShapes) ./ S_shape.vnVisit_shape(viShapes)); ylabel('t per visit (s)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
subplot 248; bar(S_shape.vpBackward_shape(viShapes)); ylabel('Backward swim prob.'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
xtickangle(hFig.Children, -30);

% Export to csv
vcAnimals = cell2mat(S_trialset.csAnimals);
vcFile_shapes_probe = subsFileExt_(vcFile_trialset, sprintf('_shapes_probe_%s.csv', vcAnimals));
[m1_,v2_,v3_,v4_] = struct_get_(S_shape, 'mrDRVS_shape', 'vnVisit_shape', 'vtVisit_shape', 'vpBackward_shape');
mrStats_shapes = [m1_;v2_;v3_;v2_./v3_;v4_]';
csStats = {'Sampling density (counts/m)', 'Sampling rate (Hz)', 'Speed (m/s)', ...
    'Escan density (counts/m)', 'Visit count', 'Visit duration (s)', 'Duration per visit', 'Freq. backward swimming'};
csvwrite_(vcFile_shapes_probe, mrStats_shapes, 'States by shapes (probe trials); [shapes, stats]');
fprintf('\tRows: %s\n', sprintf('"%s", ', S_shape.csDist_shape{:}));
fprintf('\tColumns: %s\n', sprintf('"%s", ', csStats{:}));
end %func


%--------------------------------------------------------------------------
function S_shape = pool_probe_trialset_(S_trialset, cS_trial)

cS_probe = cS_trial(S_trialset.vlProbe);
cS_shape = cell(size(cS_probe));
for i=1:numel(cS_shape)
    cS_shape{i} = nearShapes_trial_(cS_probe{i}, S_trialset.P);        
end
vS_shape = cell2mat(cS_shape)';
% fh_pool = @(vc)cell2mat({vS_shape.(vc)}');
% csName = fieldnames(vS_shape(1));
% csName = setdiff(csName, 'csDist_shape');
csName = {'mlDist_shape', 'vrD', 'vrR', 'vrS', 'vrV'};
S_shape = struct();
for i=1:numel(csName)
    eval(sprintf('S_shape.%s=cell2mat({vS_shape.%s}'');', csName{i}, csName{i}));
end
S_shape.csDist_shape = vS_shape(1).csDist_shape;
S_shape = S_shape_calc_(S_shape, S_trialset.P);
end %func


%--------------------------------------------------------------------------
function trialset_coordinates_(vcFile_trialset)
% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
% run S141106_LearningCurve_Control.m first cell
% errordlg('Not implemented yet.'); return;

[cS_trial, S_trialset] = loadShapes_trialset_(vcFile_trialset);
P = get_set_(S_trialset, 'P', load_cfg_());

% show chart
tnShapes = countShapes_trialset_(S_trialset, cS_trial);
[hFig_overview, vhImg_overview] = plot_trialset_img_(S_trialset, single(tnShapes), [0, numel(P.csShapes)]); 
set(hFig_overview, 'Name', sprintf('# Shapes: %s', vcFile_trialset));

% create a table. make it navigatable
nFiles = numel(cS_trial);
hFig_tbl = figure_new_('FigShape', ['Shape locations: ', vcFile_trialset], [0,0,.5,1]);
iTrial = 1;
hFig_tbl.UserData = makeStruct_(S_trialset, P, iTrial, tnShapes, vhImg_overview);
set0_(cS_trial);
hFig_tbl.KeyPressFcn = @(h,e)keypress_FigShape_(h,e);
plotShapes_trial_(hFig_tbl, iTrial);
% uiwait(msgbox('Right-click on the shapes and food to fill the table. Press "OK" when finished.'));
msgbox('Right-click on the shapes and food to fill the table. Close the figure when finished.');
uiwait(hFig_tbl);

% save
% if ~isvalid(hFig_tbl), msgbox('Table is closed by user, nothing is saved.'); return; end
if ~questdlg_('Save the coordinates?'), return; end
hMsgbox = msgbox('Saving... (This closes automatically)');
vcFile_mat = strrep(vcFile_trialset, '.trialset', '_trialset.mat');
save_var_(vcFile_mat, 'cS_trial', get0_('cS_trial'));
close_(hFig_tbl, hFig_overview, hMsgbox);
msgbox_(['Shape info saved to ', vcFile_mat]);
end %func


%--------------------------------------------------------------------------
% 08/12/18 JJJ: Get yes or no answer from the user
function flag = questdlg_(vcTitle, flag)
% flag: default is yes (1) or no (0)
if nargin<2, flag = 1; end
if flag
    flag = strcmpi(questdlg(vcTitle,'','Yes','No','Yes'), 'Yes');
else
    flag = strcmpi(questdlg(vcTitle,'','Yes','No','No'), 'Yes');
end
end %func


%--------------------------------------------------------------------------
function save_var_(vcFile_mat, vcName, val)
fAppend = exist_file_(vcFile_mat);
eval(sprintf('%s=val;', vcName));
if fAppend
    try
        save(vcFile_mat, vcName, '-v7.3', '-append', '-nocompression'); %faster    
    catch
        save(vcFile_mat, vcName, '-v7.3', '-append'); % backward compatible
    end
else
    try
        save(vcFile_mat, vcName, '-v7.3', '-nocompression'); %faster    
    catch
        save(vcFile_mat, vcName, '-v7.3'); % backward compatible
    end
end
end %func


%--------------------------------------------------------------------------
function plotTraj_trial_(hFig_tbl, iTrial)
S_fig = hFig_tbl.UserData;
delete_(get_(S_fig, 'hTraj'));
hTraj = plot(S_fig);
hFig_tbl.UserData = struct_add_(S_fig, hTraj);
end %func

        
%--------------------------------------------------------------------------
function plotShapes_trial_(hFig_tbl, iTrial)

S_fig = hFig_tbl.UserData;
cS_trial = get0_('cS_trial');
S_ = cS_trial{iTrial};


P1 = S_fig.P;
P1.nSkip_img = get_set_(P1, 'nSkip_img', 2);
P1.xy0 = S_.xy0 / P1.nSkip_img;
P1.pixpercm = P1.pixpercm / P1.nSkip_img;
img0 = imadjust_mask_(binned_image_(S_.img0, P1.nSkip_img));
[~,dataID_,~] = fileparts(S_.vidFname);    

% Crate axes
hAxes = get_(S_fig, 'hAxes');
if isempty(hAxes)
    hAxes = axes(hFig_tbl, 'Units', 'pixels', 'Position', [10,220,800,600]);
end

% draw figure
hImage = get_(S_fig, 'hImage');
if isempty(hImage)
    hImage = imshow(img0, 'Parent', hAxes);
    hold(hAxes, 'on');
else
    hImage.CData = img0;
end
hImage.UserData = P1;

% draw a grid
delete_(get_(S_fig, 'hGrid'));
hGrid = draw_grid_(hImage, -10:5:10);

% Title
vcTitle = [dataID_, '  [H]elp, [T]rajectory, [L/R/PgDn/PgUp]:Next/Prev, [G]oto, [E]xport ...'];
hTitle = get_(S_fig, 'hTitle');
if isempty(hTitle)
    hTitle = title_(hAxes, vcTitle);
else
    hTitle.String = vcTitle;
end

% Draw a table
hTable = get_(S_fig, 'hTable');
if isempty(hTable)
    hTable = uitable(hFig_tbl, 'Data', S_.mrPos_shape, ...
        'Position', [10 10 400 200], 'RowName', P1.csShapes, ...
        'ColumnName', {'X pos (grid)', 'Y pos (grid)', 'Angle (deg)'});
    hTable.ColumnEditable = true(1, 3);    
    hTable.CellEditCallback = @(a,b)draw_shapes_tbl_(hImage, hTable, iTrial);
else
    hTable.Data = S_.mrPos_shape;
end

% Update
delete_(get_(S_fig, 'vhShapes'));
vhShapes = draw_shapes_tbl_(hImage, hTable, iTrial);
contextmenu_(hImage, hTable);
hFig_tbl.UserData = struct_add_(S_fig, hAxes, hImage, hTable, hGrid, iTrial, vhShapes, vhShapes);
end %func


%--------------------------------------------------------------------------
function img_adj = imadjust_mask_(img, mlMask)
if nargin<2, mlMask = []; end
if isempty(mlMask)
    int_lim = quantile(single(img(img>0)), [.01, .99]);
else
    int_lim = quantile(single(img(~mlMask)), [.01, .99]);
end
% imadjust excluding the mask
img_adj = imadjust(img, (int_lim)/255, [0, 1]);
end %func


%--------------------------------------------------------------------------
function keypress_FigShape_(hFig, event)
S_fig = get(hFig, 'UserData');
nStep = 1 + key_modifier_(event, 'shift')*3;
cS_trial = get0_('cS_trial');
nTrials = numel(cS_trial);
S_trial = cS_trial{S_fig.iTrial};

switch lower(event.Key)
    case 'h'
        msgbox(...
            {'[H]elp', 
            '(Shift)+[L/R]: next trial (Shift: quick jump)', 
            '[G]oto trial', 
            '[Home]: First trial', 
            '[END]: Last trial', 
            '[E]xport coordinates to csv',
            '[T]rajectory toggle',
            '[S]ampling density', 
            '[C]opy trialset path'}, ...
                'Shortcuts');        
            
    case {'leftarrow', 'rightarrow', 'home', 'end'}
        % move to different trials and draw
        iTrial_prev = S_fig.iTrial;
        if strcmpi(event.Key, 'home')
            iTrial = 1;
        elseif strcmpi(event.Key, 'end')
            iTrial = nTrials;
        elseif strcmpi(event.Key, 'leftarrow')      
            iTrial = max(S_fig.iTrial - nStep, 1);
        elseif strcmpi(event.Key, 'rightarrow')
            iTrial = min(S_fig.iTrial + nStep, nTrials);      
        end
        if iTrial ~= iTrial_prev
            plotShapes_trial_(hFig, iTrial);
        end
        if isvalid_(get_(S_fig, 'hTraj')) % update the trajectory if turned on
            draw_traj_trial_(hFig, iTrial);
        end
        
    case 'g'
        vcTrial = inputdlg('Trial ID: ');
        vcTrial = vcTrial{1};
        if isempty(vcTrial), return; end
        vcTrial = path2DataID_(vcTrial);
        csDataID = getDataID_cS_(cS_trial);
        iTrial = find(strcmp(vcTrial, csDataID));
        if isempty(iTrial)
            msgbox(['Trial not found: ', vcTrial]);
            return; 
        end        
        plotShapes_trial_(hFig, iTrial);
        if isvalid_(get_(S_fig, 'hTraj')) % update the trajectory if turned on
            draw_traj_trial_(hFig, iTrial);
        end
        
    case 'e'
        trial2csv_(S_trial);
        
    case 't' % draw trajectory
        if isvalid_(get_(S_fig, 'hTraj'))
            delete_plot_(hFig, 'hTraj');
        else
            draw_traj_trial_(hFig, S_fig.iTrial);
        end  
        
    case 's' % Sampling density
        [S_shape, mrTXYARDVS_rs, P1] = nearShapes_trial_(S_trial, S_fig.P);        
        [vrX, vrY] = cm2pix_(mrTXYARDVS_rs(:,2:3)*100, P1);
        vrD = mrTXYARDVS_rs(:,6);
        S_trial = struct_add_(S_trial, vrX, vrY, vrD);
        hFig_grid = figure_new_('FigGrid', S_trial.vcFile_Track, [0,0,.5,.5]);
        [RGB, mrPlot] = gridMap_(S_trial, P1, 'density');
        imshow(RGB); title('Sampling density');
        
        hFig = figure_new_('',S_trial.vcFile_Track, [0 .5 .5 .5]); 
        viShapes = 1:6;
        subplot 241; bar(1./S_shape.mrDRVS_shape(1,viShapes)); ylabel('Sampling density'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
        subplot 242; bar(S_shape.mrDRVS_shape(2,viShapes)); ylabel('Sampling Rate (Hz)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
        subplot 243; bar(S_shape.mrDRVS_shape(3,viShapes)); ylabel('Speed (m/s)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
        subplot 244; bar(1./S_shape.mrDRVS_shape(4,viShapes)); ylabel('EScan density'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;        
        subplot 245; bar(S_shape.vnVisit_shape(viShapes)); ylabel('# visits'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
        subplot 246; bar(S_shape.vtVisit_shape(viShapes)); ylabel('t visit (s)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
        subplot 247; bar(S_shape.vtVisit_shape(viShapes) ./ S_shape.vnVisit_shape(viShapes)); ylabel('t per visit (s)'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;
        subplot 248; bar(S_shape.vpBackward_shape(viShapes)); ylabel('Backward swim prob.'); set(gca,'XTickLabel', S_shape.csDist_shape); grid on;

    case 'c'
        clipboard('copy', S_trial.vcFile_Track);
        msgbox(sprintf('%s copied to clipboard', S_trial.vcFile_Track));
        
    otherwise
        return;
end
end %func


%--------------------------------------------------------------------------
function [S_shape, mrTXYARDVS_rs, P1] = nearShapes_trial_(S_trial, P)
P1 = setfield(P, 'xy0', S_trial.xy0);
mrTXYARDVS_rs = resample_trial_(S_trial, P1);
[mrXY_h, vrD, vrR, vrV, vrS] = ...
    deal(mrTXYARDVS_rs(:,2:3), mrTXYARDVS_rs(:,6), mrTXYARDVS_rs(:,5), mrTXYARDVS_rs(:,7), mrTXYARDVS_rs(:,8));
% vrV = hypot(diff3_(mrXY_h(:,1)), diff3_(mrXY_h(:,2))) * get_set_(P, 'sRateHz_resample', 100); % meter per sec

cm_per_grid = get_set_(P, 'cm_per_grid', 5);
mrPos_shape_meter = S_trial.mrPos_shape;
mrPos_shape_meter(:,1:2) = mrPos_shape_meter(:,1:2) * cm_per_grid / 100;
nShapes = size(mrPos_shape_meter,1);
mlDist_shape = false(size(mrXY_h,1), nShapes);
dist_cut = get_set_(P, 'dist_cm_shapes', 3) / 100;

for iShape = 1:nShapes
    vcShape = strtok(P.csShapes{iShape}, ' ');
    xya_ = mrPos_shape_meter(iShape,:);
    len_ = P.vrShapes(iShape)/100; % in meter
    [mrXY_poly_, fCircle] = get_polygon_(vcShape, xya_(1:2), len_, xya_(3));
    if fCircle
        vrD_ = hypot(xya_(1)-mrXY_h(:,1), xya_(2)-mrXY_h(:,2)) - len_/2; 
    else %polygon
        vrD_ = nearest_perimeter_(mrXY_poly_, mrXY_h); % convert to meter
    end
    mlDist_shape(:,iShape) = vrD_ <= dist_cut;
end

% distance to the wall
r_wall = get_set_(P, 'diameter_cm_wall', 150) / 2 / 100;
dist_wall = get_set_(P, 'dist_cm_wall', 15) / 100; % 15 cm from the wall
vl_wall = hypot(mrXY_h(:,1), mrXY_h(:,2)) >= (r_wall - dist_wall);
mlDist_shape = [mlDist_shape, vl_wall, ~vl_wall];
csDist_shape = [P.csShapes, 'Wall', 'Not Wall'];

S_shape = makeStruct_(mlDist_shape, csDist_shape, vrD, vrR, vrV, vrS);
S_shape = S_shape_calc_(S_shape, P);
end %func


%--------------------------------------------------------------------------
function S_shape = S_shape_calc_(S_shape, P)
% output 
% ------
% mrDRVS_shape
% vnVisit_shape
% vtVisit_shape
% vpBackward_shape
[vrD, vrR, vrV, vrS, mlDist_shape] = struct_get_(S_shape, 'vrD', 'vrR', 'vrV', 'vrS', 'mlDist_shape');

% sampling density by shapes
sRateHz_rs = get_set_(P, 'sRateHz_resample', 100);
mrDRVS_shape = region_median_([vrD, vrR, abs(vrV), vrS], mlDist_shape, @nanmedian); %@nanmean
[~, vnVisit_shape] = findup_ml_(mlDist_shape, sRateHz_rs);
% vnVisit_shape = sum(diff(mlDist_shape)>0); % clean up transitions
vtVisit_shape = sum(mlDist_shape) / sRateHz_rs;
vpBackward_shape = region_median_(sign(vrV)<0, mlDist_shape, @nanmean);

S_shape = struct_add_(S_shape, mrDRVS_shape, vnVisit_shape, vtVisit_shape, vpBackward_shape);
end %func


%--------------------------------------------------------------------------
function [cvi, vn] = findup_ml_(ml, nRefrac)
cvi = cell(1, size(ml,2));
vn = zeros(1, size(ml,2));
for iCol=1:size(ml,2)
    vi_ = find(diff(ml(:,iCol))>0);
    vn_ = diff(vi_);
    viiKill_ = find(vn_<nRefrac);
    if ~isempty(viiKill_)
        vi_(viiKill_ + 1) = []; % remove 
    end
    cvi{iCol} = vi_;
    vn(iCol) = numel(vi_);
end
% vn = cell2mat(cellfun(@(x)diff(x), cvi', 'UniformOutput', 0));
end %func


%--------------------------------------------------------------------------
function mrMed = region_median_(mr, ml, fh)
if nargin<3, fh = @median; end
mrMed = nan(size(mr,2), size(ml,2));
for iCol = 1:size(ml,2)
    vi_ = find(ml(:,iCol));
    if ~isempty(vi_)
        mrMed(:,iCol) = fh(mr(vi_,:));
    end
end
end %func


%--------------------------------------------------------------------------
function S_fig = delete_plot_(hFig, vcTag)
S_fig = hFig.UserData;
if isempty(S_fig), return; end
delete_(get_(S_fig, vcTag));
S_fig.(vcTag) = [];
hFig.UserData = S_fig;
end %func


%--------------------------------------------------------------------------
function [S_fig, hPlot] = draw_traj_trial_(hFig, iTrial)
S_fig = hFig.UserData;
cS_trial = get0_('cS_trial');
S_trial = cS_trial{iTrial};
P = get_set_(S_fig, 'P', load_cfg_());
nSkip_img = get_set_(P, 'nSkip_img', 2);

[S_shape, mrTXYARDVS_rs, P1] = nearShapes_trial_(S_trial, P);
mrXY_pix = cm2pix_(mrTXYARDVS_rs(:,2:3)*100, P1);

try
%     [X,Y] = deal(S_trial.XC(:,2)/nSkip_img, S_trial.YC(:,2)/nSkip_img);
    [X,Y] = deal(mrXY_pix(:,1)/nSkip_img, mrXY_pix(:,2)/nSkip_img);
%     vl = ~S_shape.mlDist_shape(:,end); % location query
%     [X(vl),Y(vl)] = deal(nan(sum(vl),1));
    hPlot = get_(S_fig, 'hTraj');
    if isvalid_(hPlot)
        [hPlot.XData, hPlot.YData] = deal(X, Y);    
    else
        S_fig.hTraj = plot(S_fig.hAxes, X, Y, 'b');
        hFig.UserData = S_fig;
    end
catch
    ; % pass
end
end %func


%--------------------------------------------------------------------------
% Check for shift, alt, ctrl press
function flag = key_modifier_(event, vcKey)
try
    flag = any(strcmpi(event.Modifier, vcKey));
catch
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
% Count number of shapes input
function tnShapes = countShapes_trialset_(S_trialset, cS_trial)
[viImg, tiImg] = struct_get_(S_trialset, 'viImg', 'tiImg');
tnShapes = zeros(size(tiImg), 'uint8');
for iTrial = 1:numel(cS_trial)    
    S_ = cS_trial{iTrial};
    mrPos_shape = get_(S_, 'mrPos_shape');
    if ~isempty(mrPos_shape)
        tnShapes(viImg(iTrial)) = sum(~any(isnan(mrPos_shape), 2));
    end
end %for
end %func


%--------------------------------------------------------------------------
function contextmenu_(hImg, tbl)
c = uicontextmenu;
hImg.UIContextMenu = c;
P1 = hImg.UserData;
% Create child menu items for the uicontextmenu
csShapes = P1.csShapes;
for iShape=1:numel(csShapes)
    uimenu(c, 'Label', csShapes{iShape}, 'Callback',@setTable_);
end
uimenu(c, 'Label', '--------');
uimenu(c, 'Label', 'Rotate CW', 'Callback',@setTable_);
uimenu(c, 'Label', 'Rotate CCW', 'Callback',@setTable_);
uimenu(c, 'Label', 'Delete', 'Callback',@setTable_);

    function setTable_(source,callbackdata)   
        xy = get(hImg.Parent, 'CurrentPoint');
        xy_cm = pix2cm_(xy(1,1:2), P1); % scale factor
        xy_grid = round(xy_cm / P1.cm_per_grid); 
        iRow_nearest = findNearest_grid_(xy_grid, tbl.Data, 1);
        switch lower(source.Label)
            case {'rotate cw', 'rotate ccw'}                
                if isempty(iRow_nearest), return; end
                dAng = ifeq_(strcmpi(source.Label, 'rotate cw'), 90, -90);
                tbl.Data(iRow_nearest,3) = mod(tbl.Data(iRow_nearest,3)+dAng,360);
            case 'delete'
                if isempty(iRow_nearest), return; end
                tbl.Data(iRow_nearest,:) = nan; %delete
            otherwise
                if ~isempty(iRow_nearest)
                    tbl.Data(iRow_nearest,:) = nan; %delete
                end
                iRow = find(strcmp(tbl.RowName, source.Label));
                tbl.Data(iRow,:) = [xy_grid(:)', 0];
        end %switch
        draw_shapes_tbl_(hImg, tbl);
    end
end %func


%--------------------------------------------------------------------------
function iRow_nearest = findNearest_grid_(xy_grid, mrGrid, d_max);
d = pdist2(xy_grid(:)', mrGrid(:,1:2));
iRow_nearest = find(d<=d_max, 1, 'first');
end %func


%--------------------------------------------------------------------------
function h = draw_grid_(hImg, viGrid)
P1 = hImg.UserData;
[xx_cm, yy_cm] = meshgrid(viGrid);
mrXY_pix = cm2pix_([xx_cm(:), yy_cm(:)] * P1.cm_per_grid, P1);
h = plot(hImg.Parent, mrXY_pix(:,1), mrXY_pix(:,2), 'r.');
end %func


%--------------------------------------------------------------------------
function vhPlot = draw_shapes_tbl_(hImg, tbl, iTrial)

hFig = hImg.Parent.Parent;
S_fig = hFig.UserData;
if nargin<3, iTrial = S_fig.iTrial; end
P1 = hImg.UserData;
delete_(get_(P1, 'vhPlot'));
mrXY = tbl.Data;
mrXY(:,1:2) = mrXY(:,1:2) * P1.cm_per_grid;
nShapes = size(mrXY,1);
vhPlot = zeros(nShapes, 1);
for iShape = 1:nShapes
    vcShape_ = strtok(tbl.RowName{iShape}, ' ');
    vhPlot(iShape) = draw_shapes_img_(hImg, mrXY(iShape,:), P1.vrShapes(iShape), vcShape_);
end
P1.vhPlot = vhPlot;
set(hImg, 'UserData', P1);

% Save table data to fig userdata
cS_trial = get0_('cS_trial');
cS_trial{iTrial}.mrPos_shape = tbl.Data;
set0_(cS_trial);
S_fig.vhShapes = vhPlot;
hFig.UserData = S_fig;
end %func


%--------------------------------------------------------------------------
function flag = isvalid_(h)
if isempty(h), flag = 0; return ;end
try
    flag = isvalid(h);
catch
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function h = draw_shapes_img_(hImg, xya, dimm, vcShape)
% xya: xy0 and angle (cm and deg)
h = nan; 
if any(isnan(dimm)), return; end
xy_ = xya(1:2);
if any(isnan(xy_)), return; end

P1 = hImg.UserData;
if numel(xya)==3
    ang = xya(3); 
else
    ang = 0;
end

mrXY_cm = get_polygon_(vcShape, xy_, dimm, ang);
mrXY_pix = cm2pix_(mrXY_cm, P1);
h = plot(hImg.Parent, mrXY_pix(:,1), mrXY_pix(:,2), 'g-', 'LineWidth', 1);
end %func


%--------------------------------------------------------------------------
function [mrXY_cm, fCircle] = get_polygon_(vcShape, xy_, dimm, ang)
fCircle = 0;
switch upper(vcShape)
    case 'TRIANGLE' % length is given
        r_ = dimm(1)/sqrt(3);
        vrA_ = [0, 120, 240, 0];
    case {'CIRCLE', 'FOOD'} % diameter is given
        r_ = dimm(1)/2;
        vrA_ = [0:9:360]; 
        fCircle = 1;
    case {'SQUARE', 'RECT', 'RECTANGLE'} % length is given
        r_ = dimm(1);
        vrA_ = [45:90:360+45];
    otherwise, error(['draw_shapes_img_: invalid shape: ', vcShape]);
end %switch
mrXY_cm = bsxfun(@plus, xy_(:)', rotate_line_(vrA_ + ang, r_));
end %func


%--------------------------------------------------------------------------
function xy_cm = pix2cm_(xy_pix, P1)
xy_cm = bsxfun(@minus, xy_pix, P1.xy0(:)') / P1.pixpercm;
xy_cm(:,2) = -xy_cm(:,2); % image coordinate to xy coordinate
xy_cm = rotatexy_(xy_cm, -P1.angXaxis);
end %func


%--------------------------------------------------------------------------
function varargout = cm2pix_(xy_cm, P1)
xy_pix = xy_cm * P1.pixpercm;
xy_pix(:,2) = -xy_pix(:,2); % change y axis
xy_pix = bsxfun(@plus, rotatexy_(xy_pix, -P1.angXaxis), P1.xy0(:)');
if nargout==1
    varargout{1} = xy_pix;
else
    [varargout{1}, varargout{2}] = deal(xy_pix(:,1), xy_pix(:,2));
end
end %func


%--------------------------------------------------------------------------
function [ xy_rot ] = rotatexy_( xy, ang)
%ROTATEXY rotate a vector with respect to the origin, ang in degree
% xy = xy(:);
CosA = cos(deg2rad(ang));
SinA = sin(deg2rad(ang));
M = [CosA, -SinA; SinA, CosA];
xy_rot = (M * xy')';
end


%--------------------------------------------------------------------------
% rotate a line and project. rotate from North
function xy = rotate_line_(vrA_deg, r)
if nargin<2, r=1; end
vrA_ = pi/2 - vrA_deg(:)/180*pi;
xy = r * [cos(vrA_), sin(vrA_)];
end %func


%--------------------------------------------------------------------------
function img1 = binned_image_(img, nSkip, fFast)
% fFast: set to 0 to do averaging (higher image quality)
if nargin<3, fFast = 1; end
if ndims(img)==3, img = img(:,:,1); end
if fFast
    img1 = img(1:nSkip:end, 1:nSkip:end); % faster
else
    dimm1 = floor(size(img)/nSkip);
    viY = (0:dimm1(1)-1) * nSkip;
    viX = (0:dimm1(2)-1) * nSkip;
    img1 = zeros(dimm1, 'single');
    for ix = 1:nSkip
        for iy = 1:nSkip
            img1 = img1 + single(img(viY+iy, viX+ix));
        end
    end
    img1 = img1 / (nSkip*nSkip);
    if isa(img, 'uint8'), img1 = uint8(img1); end
end
end %func


%--------------------------------------------------------------------------
function [cS_trial, S_trialset, mrPath, mrDur] = loadShapes_trialset_(vcFile_trialset)

[mrPath, mrDur, S_trialset, cS_trial] = trialset_learningcurve_(vcFile_trialset);
nSkip_img = get_set_(S_trialset.P, 'nSkip_img', 2);
% if nargout>=3
%     trImg0 = cellfun(@(x)imadjust(binned_image_(x.img0, nSkip_img)), cS_trial, 'UniformOutput', 0);
%     trImg0 = cat(3, trImg0{:});
% end

% default shape table
csShapes = get_set_(S_trialset, 'csShapes', {'Triangle Lg', 'Triangle Sm', 'Square Lg', 'Square Sm', 'Circle Lg', 'Circle Sm', 'Food'});
csShapes = csShapes(:);
nShapes = numel(csShapes);
mrData0 = [nan(nShapes, 2), zeros(nShapes,1)];

% load prev result
vcFile_mat = strrep(vcFile_trialset, '.trialset', '_trialset.mat');
[cTable_data, cS_trial_prev] = load_mat_(vcFile_mat, 'cTable_data', 'cS_trial');

% fill in mrPos_shape
csDataID = getDataID_cS_(cS_trial);
csDataID_prev = getDataID_cS_(cS_trial_prev);
for iFile = 1:numel(cS_trial)
    S_ = cS_trial{iFile};
    if isfield(S_, 'mrPos_shape'), continue; end
    iPrev = find(strcmp(csDataID{iFile}, csDataID_prev));
    mrData_ = mrData0;
    if ~isempty(iPrev)     
        mrData_prev = cS_trial_prev{iPrev}.mrPos_shape;
        nCol = min(nShapes, size(mrData_prev,1));
        mrData_(1:nCol,:) = mrData_prev(1:nCol,:);
    end
    S_.mrPos_shape = mrData_;    
    cS_trial{iFile} = S_;
end
end %func


%--------------------------------------------------------------------------
function csDataID = getDataID_cS_(cS)
csDataID = cell(size(cS));
for i=1:numel(cS)
    [~,csDataID{i},~] = fileparts(cS{i}.vidFname);
end
end %func


%--------------------------------------------------------------------------
function dataID = path2DataID_(vc)
if iscell(vc), vc = vc{1}; end
[~, dataID, ~] = fileparts(vc);
dataID = strrep(dataID, '_Track', '');
end %func


%--------------------------------------------------------------------------
function h = msgbox_(vcMsg, fEcho)
if nargin<2, fEcho = 1; end
h = msgbox(vcMsg);
if fEcho, disp(vcMsg); end
end %func


%--------------------------------------------------------------------------
% 7/26/2018 JJJ: save mat file
% function save_mat_(varargin)
% vcFile = varargin{1};
% for i=2:nargin
%     eval('%s=varargin{%d};', inputname(i));
% end
% if exist_file_(vcFile)    
%     save(vcFile, varargin{2:end}, '-append');
% else
%     save(vcFile, varargin{2:end});
% end
% end %func


%--------------------------------------------------------------------------
function varargout = load_mat_(varargin)
if nargin<1, return; end
vcFile_mat = varargin{1};
varargout = cell(1, nargout());
if ~exist_file_(vcFile_mat), return; end
if nargin==1, S = load(vcFile_mat); return; end

S = load(vcFile_mat, varargin{2:end});
for iArg = 1:nargout()
    try
        varargout{iArg} = getfield(S, varargin{iArg+1});
    catch
        ;    
    end
end %for
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
S_cfg.vcFile_settings = get_set_(S_cfg, 'vcFile_settings', 'settings_vistrack.m');
S_cfg.pixpercm = get_set_(S_cfg, 'pixpercm', 7.238);
S_cfg.angXaxis = get_set_(S_cfg, 'angXaxis', -0.946);
end %func


%--------------------------------------------------------------------------
function trialset_exportcsv_(vcFile_trialset)
csMsg = {'Exporting the trialset to csv files...(this will close when done)', 'It can take up to several minutes'};
h = msgbox(csMsg, 'modal');
% S_trialset = load_trialset_(vcFile_trialset);
[cS_trial, S_trialset] = loadShapes_trialset_(vcFile_trialset);
% csFiles_track = S_trialset.csFiles_Track;
% csFiles_failed = {};
% for iFile = 1:numel(csFiles_track)
for iFile = 1:numel(cS_trial)
    S_ = cS_trial{iFile};
    if isempty(S_), continue; end
    try
        [~,~,vcMsg,csFormat] = trial2csv_(S_, S_trialset.P);
        fprintf('%s\n', vcMsg);
    catch
        disp(lasterr());
    end
end %for
disp_cs_(csFormat);
close_(h);
end %func


%--------------------------------------------------------------------------
function trial_gridmap_(vcFile_Track)
S_trial = load_(vcFile_Track);
P = load_settings_(S_trial);
% LOADSETTINGS;
h = msgbox('Calculating... (this will close automatically)');

S_ = importTrial(S_trial, P.pixpercm, P.angXaxis);
[RGB, mrPlot] = gridMap_(S_, P, 'time');
% [mnVisit1, mnVisit] = calcVisitCount(S_, S_.img0);
% dataID = S_trial
figure_new_('', S_trial.vidFname); imshow(RGB);
title('Time spent');
try close(h); catch, end;   
end %func


%--------------------------------------------------------------------------
function trial_timemap_(S_trial)
P = load_settings_(S_trial);
%track head
h = msgbox('Calculating... (This closes automatically)');
[VISITCNT, TIMECNT] = calcVisitDensity(S_trial.img0, S_trial.TC, S_trial.XC(:,2), S_trial.YC(:,2), P.TRAJ_NFILT);

% trialID = trial_id_(handles);
img0_adj = imadjust(S_trial.img0);
hFig = figure_new_('', S_trial.vidFname); 
imshow(rgbmix_(img0_adj, TIMECNT));
resize_figure_(hFig, [0,0,.5,1]);
title('Time map');
close_(h);
end %func


%--------------------------------------------------------------------------
function [dataID, fishID, iSession, iTrial, fProbe] = trial_id_(S_trial)
[~,dataID,~] = fileparts(S_trial.vidFname);
fishID = dataID(4);
iSession = str2num(dataID(2:3));
iTrial = str2num(dataID(5));
fProbe = numel(dataID) > 5;
end %func


%--------------------------------------------------------------------------
function [RGB, mrPlot] = gridMap_(vsTrial, P, vcMode, lim, mlMask)
% vcMode: {'time', 'visit', 'time/visit', 'density'}

if nargin < 2, P = []; end
if nargin < 3, vcMode = 'time'; end % visit, time, time/visit
if nargin < 4, lim = []; end
if nargin < 5, mlMask = []; end

nGrid_map = get_set_(P, 'nGrid_map', 20);
nTime_map = get_set_(P, 'nTime_map', 25);
angXaxis = get_set_(P, 'angXaxis', -1.1590); %deg

if iscell(vsTrial), vsTrial = cell2mat(vsTrial); end % make it an array

%background image processing
xy0 = vsTrial(1).xy0;
img0 = vsTrial(1).img0(:,:,1);    
% mlMask = getImageMask(img0, [0 60], 'CENTRE');
img0 = imrotate(imadjust(img0), -angXaxis, 'nearest', 'crop');

%rotate vrX, vrY, and images
vrX = poolVecFromStruct(vsTrial, 'vrX'); % in meters
vrY = poolVecFromStruct(vsTrial, 'vrY'); %  in meters
try vrD = poolVecFromStruct(vsTrial, 'vrD'); catch, vrD = []; end
rotMat = rotz(-angXaxis);    rotMat = rotMat(1:2, 1:2);
mrXY = [vrX(:) - xy0(1), vrY(:) - xy0(2)] * rotMat;
vrX = mrXY(:,1) + xy0(1);
vrY = mrXY(:,2) + xy0(2);

viX = ceil(vrX/nGrid_map);
viY = ceil(vrY/nGrid_map);
[h, w] = size(img0);
h = h / nGrid_map;
w = w / nGrid_map;
[mrDensity, mnVisit, mnTime] = deal(zeros(h,w));
for iy=1:h
    vlY = (viY == iy);
    for ix=1:w
        viVisit = find(vlY & (viX == ix));  
        if isempty(viVisit), continue; end
        mnTime(iy,ix) = numel(viVisit);
        nRepeats = sum(diff(viVisit) < nTime_map); % remove repeated counts
        mnVisit(iy,ix) = numel(viVisit) - nRepeats;      
        if ~isempty(vrD)
            mrDensity(iy,ix) = 1 ./ mean(vrD(viVisit));
            fprintf('.');
        end
    end
end

mrTperV = mnTime ./ mnVisit;

switch lower(vcMode)
    case 'time'
        mrPlot = mnTime;
    case 'visit'
        mrPlot = mnVisit;
    case 'time/visit'
        mrPlot = mrTperV;
    case {'density', 'samplingdensity'}
        mrPlot = mrDensity;
end

mnPlot_ = imresize(mrPlot, nGrid_map, 'nearest');
if isempty(lim), lim = [min(mnPlot_(:)) max(mnPlot_(:))]; end
mrVisit = uint8((mnPlot_ - lim(1)) / diff(lim) * 255);
RGB = rgbmix_(img0, mrVisit, mlMask);
if nargout==0
    figure; imshow(RGB); title(sprintf('%s, clim=[%f, %f]', vcMode, lim(1), lim(2)));
end
end %func


%--------------------------------------------------------------------------
function img = rgbmix_(img_bk, img, MASK, mixRatio)

if nargin<3, MASK = []; end
if nargin<4, mixRatio = []; end

if isempty(mixRatio), mixRatio = .25; end

if numel(size(img_bk)) == 2 %gray scale
    if ~isa(img_bk, 'uint8')
        img_bk = uint8(img_bk/max(img_bk(:)));
    end
    img_bk = imgray2rgb(img_bk, [0 255], 'gray');
end
if numel(size(img)) == 2 %gray scale
    if ~isempty(MASK), img(~MASK) = 0; end % clear non masked area (black)
    if ~isa(img, 'uint8')
        img = uint8(img/max(img(:))*255);
    end
    img = imgray2rgb(img, [0 255], 'jet');
end

for iColor = 1:3
    mr1_ = single(img(:,:,iColor));
    mr0_ = single(img_bk(:,:,iColor));
    mr_ = mr1_*mixRatio + mr0_*(1-mixRatio);            
    if isempty(MASK)
        img(:,:,iColor) = uint8(mr_);
    else
        mr0_(MASK) = mr_(MASK);
        img(:,:,iColor) = uint8(mr0_);
    end
end
end %func


%--------------------------------------------------------------------------
function handles = trial_fixsync1_(handles, fAsk)
% Load video file from handle
if nargin<2, fAsk = 1; end

% Load video
h=msgbox('Loading... (this will close automatically)');
[vidobj, vcFile_vid] = load_vid_handle_(handles);
if isempty(vidobj)
    fprintf(2, 'Video file does not exist: %s\n', handles.vidFname);
    close_(h);
    return;
end
P = load_settings_(handles);

% load video, load LED until end of the video
try
    nFrames_load = handles.FLIM(2);
catch
    nFrames_load = vidobj.NumberOfFrames;
end
[vrLed_cam, viT_cam] = loadLed_vid_(vidobj, [], nFrames_load);
[viT_cam, viiFilled_led] = fill_pulses_(viT_cam);
close_(h);
% figure; plot(vrLed); hold on; plot(viT_cam, vrLed(viT_cam), 'o'); 


% get ADC timestamp
vrT_adc = getSync_adc_(handles);
nBlinks = min(numel(viT_cam), numel(vrT_adc));
[viT_cam, vrT_adc] = deal(viT_cam(1:nBlinks), vrT_adc(1:nBlinks));

% Compare errors
vtLed_cam = interp1(viT_cam, vrT_adc, (1:numel(vrLed_cam)), 'linear', 'extrap');
[vrX, vrY, TC] = deal(handles.XC(:,2), handles.YC(:,2), handles.TC(:));
vrTC_new = interp1(viT_cam, vrT_adc, (handles.FLIM(1):handles.FLIM(2))', 'linear', 'extrap');
vrT_err = TC - vrTC_new;
vrV = sqrt((vrX(3:end)-vrX(1:end-2)).^2 + (vrY(3:end)-vrY(1:end-2)).^2) / P.pixpercm / 100;
vrV_prev = vrV ./ (TC(3:end) - TC(1:end-2));
vrV_new = vrV ./ (vrTC_new(3:end) - vrTC_new(1:end-2));

% plot
hFig = figure_new_('', vcFile_vid); 
ax(1) = subplot(3,1,1);
plot(vtLed_cam, vrLed_cam); grid on; hold on;
plot(vtLed_cam(viT_cam), vrLed_cam(viT_cam), 'ro');
ylabel('LED');
title(sprintf('FPS: %0.3f Hz', handles.FPS));

ax(2) = subplot(3,1,2);
plot(vrTC_new, vrT_err, 'r.'); grid on;
title(sprintf('Sync error SD: %0.3fs', std(vrT_err))); 

ax(3) = subplot(3,1,3); hold on;
plot(vrTC_new(2:end-1), vrV_prev, 'ro-'); 
plot(vrTC_new(2:end-1), vrV_new, 'go-'); grid on; 
ylabel('Speed (m/s)'); xlabel('Time (s)');
linkaxes(ax,'x');
xlim(vrTC_new([1, end]));
title(sprintf('Ave speed: %0.3f(old), %0.3f(new) m/s', mean(vrV_prev), mean(vrV_new)));

if fAsk
    vcAns = questdlg('Save time sync?', vcFile_vid, ifeq_(std(vrT_err) > .01, 'Yes', 'No'));
    fSave = strcmpi(vcAns, 'Yes');
else
    fSave = 1;
end
if fSave % save to file
    handles.TC = vrTC_new;
    handles.FPS = diff(handles.FLIM([1,end])) / diff(handles.TC([1,end]));
    trial_save_(handles);
end
end %func


%--------------------------------------------------------------------------
function hPlot = plot_vline_(hAx, vrX, ylim1, lineStyle)
if nargin<4, lineStyle = []; end
mrX = repmat(vrX(:)', [3,1]);
mrY = nan(size(mrX));
mrY(1,:) = ylim1(1);
mrY(2,:) = ylim1(2);
if isempty(lineStyle)
    hPlot = plot(hAx, mrX(:), mrY(:));
else
    hPlot = plot(hAx, mrX(:), mrY(:), lineStyle);
end
end %func


%--------------------------------------------------------------------------
function keypress_FigSync_(hFig, event)
S_fig = get(hFig, 'UserData');
if key_modifier_(event, 'shift')
    nStep = 10;
elseif key_modifier_(event, 'control')
    nStep = 100; 
else
    nStep = 1;
end
nFrames = size(S_fig.mov, 3);
switch lower(event.Key)
    case 'h'
        msgbox(...
            {'[H]elp', 
            '(Shift/Ctrl)+[L/R]: Next Frame (Shift:10x, Ctrl:100x)', 
            '[PgDn/PgUp]: Next/Prev Event Marker'            
            '[G]oto trial', 
            '[Home]: First trial', 
            '[END]: Last trial'
            }, ...
                'Shortcuts');     
        return;
        
    case {'leftarrow', 'rightarrow', 'home', 'end', 'pagedown', 'pageup'}
        % move to different trials and draw
        iFrame_prev = S_fig.iFrame;
        if strcmpi(event.Key, 'home')
            iFrame = 1;
        elseif strcmpi(event.Key, 'end')
            iFrame = nFrames;
        elseif strcmpi(event.Key, 'leftarrow')            
            iFrame = S_fig.iFrame - nStep;
        elseif strcmpi(event.Key, 'rightarrow')
            iFrame = S_fig.iFrame + nStep;
        elseif strcmpi(event.Key, 'pageup')
            iFrame = find_event_sync_(S_fig, 0);
        elseif strcmpi(event.Key, 'pagedown')            
            iFrame = find_event_sync_(S_fig, 1);
        end
        iFrame = setlim_(iFrame, [1, nFrames]);
        if iFrame ~= iFrame_prev
            refresh_FigSync_(hFig, iFrame);
        end
        
    case 'g'
        vcFrame = inputdlg('Frame#: ');
        if isempty(vcFrame), return; end
        iFrame = str2num(vcFrame);
        if isempty(iFrame) || isnan(iFrame)
            msgbox(['Invalid Frame#: ', vcFrame]);
            return; 
        end        
        refresh_FigSync_(hFig, iFrame);
        
    otherwise
        return;
end %switch
end %func


%--------------------------------------------------------------------------
function iFrame = find_event_sync_(S_fig, fForward)
if nargin<2, fForward = 1; end

% find event 
iFrame_now = S_fig.iFrame;
viText_cam = adc2cam_sync_([], S_fig.vtText);
if fForward
    iText = find(viText_cam > iFrame_now, 1, 'first');
else
    iText = find(viText_cam < iFrame_now, 1, 'last');
end
iFrame = ifeq_(isempty(iText), iFrame_now, viText_cam(iText));
if iFrame<1, iFrame = iFrame_now; end
end %func


%--------------------------------------------------------------------------
function [vtText, csText] = getText_adc_(handles, P)
if nargin<2, P=[]; end
if isempty(P), P = load_settings_(handles); end

ADCTS = get_(handles, 'ADCTS');
if isempty(ADCTS), vtText = []; return; end
ADC_CH_TEXT = get_set_(P, 'ADC_CH_TEXT', 30);
S_text = getfield(ADCTS, sprintf('%s_Ch%d', getSpike2Prefix_(ADCTS), ADC_CH_TEXT));
[vtText, vcText_] = struct_get_(S_text, 'times', 'text');
csText = cellstr(vcText_);
end %func


%--------------------------------------------------------------------------
function [vtEodr, vrEodr] = getEodr_adc_(handles, P)
if nargin<2, P=[]; end
if isempty(P), P = load_settings_(handles); end
[vtEodr, vrEodr] = deal([]);
ADCTS = get_(handles, 'ADCTS');
if isempty(ADCTS), return; end
ADC_CH_EOD = get_set_(P, 'ADC_CH_EOD', 10);
S_eod = getfield(ADCTS, sprintf('%s_Ch%d', getSpike2Prefix_(ADCTS), ADC_CH_EOD));
vtEod = S_eod.times;
vrEodr = 2 ./ (vtEod(3:end) - vtEod(1:end-2));
vtEodr = vtEod(2:end-1);
end %func


%--------------------------------------------------------------------------
function [vrLed, viT_cam] = loadLed_vid_(vidobj, xyLed, nFrames)
if nargin<2, xyLed = []; end
if nargin<3, nFrames = []; end
nStep = 300;
nParfor = 4;

t1=tic;

% Find LED
if isempty(nFrames), nFrames = vidobj.NumberOfFrames; end
flim = [1,min(nStep,nFrames)];
mov_ = vid_read(vidobj, flim(1):flim(2));
if isempty(xyLed), xyLed = findLed_mov_(mov_); end
vrLed = mov2led_(mov_, xyLed);
if flim(2) == nFrames, return; end

% Load rest of the movie
viFrame_start = (1:nStep:nFrames)';
cvrLed = cell(size(viFrame_start));
cvrLed{1} = vrLed;
try
    parfor (i = 2:numel(cvrLed), nParfor)
        flim_ = viFrame_start(i) + [0, nStep-1];
        flim_(2) = min(flim_(2), nFrames);
        cvrLed{i} = mov2led_(vid_read(vidobj, flim_(1):flim_(2)), xyLed);
    end
catch
    for iFrame1 = 2:numel(cvrLed)
        flim_ = viFrame_start(i) + [0, nStep-1];
        flim_(2) = min(flim_(2), nFrames);
        cvrLed{i} = mov2led_(vid_read(vidobj, flim_(1):flim_(2)), xyLed);
    end
end
vrLed = cell2mat(cvrLed);
if nargout>=2
    thresh_led = (max(vrLed) + median(vrLed))/2;
    viT_cam = find(diff(vrLed > thresh_led)>0) + 1;
end
fprintf('LED loading took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
% Remove pulses
function [viT_new, viRemoved] = remove_pulses_(viT)
% remove pulses out of the range
tol = .01; % allow tolerence
int_med = median(diff(viT));
int_lim = int_med * [1-tol, 1+tol];
viInt2 = viT(3:end) - viT(1:end-2);
viRemoved = find(viInt2 >= int_lim(1) & viInt2 <= int_lim(2))+1;
viT_new = viT;
if ~isempty(viRemoved)
    viT_new(viRemoved) = [];
    fprintf('Removed %d ADC pulses\n', numel(viRemoved));
end
end %func


%--------------------------------------------------------------------------
% Fill missing LED pulses
function [viT_new, viT_missing] = fill_pulses_(viT_cam)

% vlPulse = false(1, numel(viT_cam));
% vlPulse(viT_cam) = 1;
viT_ = [0; viT_cam(:)];
vrTd = diff(viT_);
vnInsert_missing = round(vrTd / median(vrTd)) - 1;
viMissing = find(vnInsert_missing > 0);
if isempty(viMissing)
    viT_new = viT_cam;
    viT_missing = [];
else
    cviT_missing = cell(1, numel(viMissing));
    for iMissing1 = 1:numel(viMissing)
        iMissing = viMissing(iMissing1);
        n_ = vnInsert_missing(iMissing);
        vi_ = linspace(viT_(iMissing), viT_(iMissing+1), n_+2);
        cviT_missing{iMissing1} = vi_(2:end-1);
    end
    viT_missing = round(cell2mat(cviT_missing));
    viT_new = sort([viT_cam(:); viT_missing(:)]); 
end
if numel(viT_missing)>0
    fprintf('%d pulses inserted (before: %d, after: %d)\n', numel(viT_missing), numel(viT_cam), numel(viT_new));
end
if nargout==0    
    figure; hold on; grid on;
    plot(viT_cam, ones(size(viT_cam)), 'bo');
    plot(viT_missing, ones(size(viT_missing)), 'ro');
end
end %func


%--------------------------------------------------------------------------
function vrLed = mov2led_(mov, xyLed)
vrLed = squeeze(mean(mean(mov(xyLed(2)+[-1:1], xyLed(1)+[-1:1], :),1),2));
end %func


%--------------------------------------------------------------------------
function vrT_adc = getSync_adc_(handles, P)
if nargin<2, P=[]; end
if isempty(P), P = load_settings_(handles); end

ADCTS = get_(handles, 'ADCTS');
if isempty(ADCTS), vrT_adc = []; return; end
S_adc = getfield(ADCTS, sprintf('%s_Ch%d', getSpike2Prefix_(ADCTS), P.ADC_CH_TCAM));
vrT_adc = get_(S_adc, 'times');
end %func


%--------------------------------------------------------------------------
function [vidobj, vcFile_vid] = load_vid_handle_(handles);
vidobj = [];
vcFile_vid = handles.vidFname;
if ~exist_file_(vcFile_vid)
    try
        vcFile_Track = get_(handles.editResultFile, 'String');
    catch
        vcFile_Track = get_(handles, 'vcFile_Track');
    end
    vcFile_vid_ = subsDir_(vcFile_vid, vcFile_Track);
    if ~exist_file_(vcFile_vid_)
        return;
    else
        vcFile_vid = vcFile_vid_;
    end
end
vidobj = get_(handles, 'vidobj');
if isempty(vidobj)
    vidobj = VideoReader_(vcFile_vid);
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_new = subsDir_(vcFile, vcDir_new)
% vcFile_new = subsDir_(vcFile, vcFile_copyfrom)
% vcFile_new = subsDir_(vcFile, vcDir_copyfrom)

% Substitute dir
if isempty(vcDir_new), vcFile_new = vcFile; return; end

[vcDir_new,~,~] = fileparts(vcDir_new); % extrect directory part. danger if the last filesep() doesn't exist
[vcDir, vcFile, vcExt] = fileparts(vcFile);
vcFile_new = fullfile(vcDir_new, [vcFile, vcExt]);
end % func


%--------------------------------------------------------------------------
function xyLed = findLed_mov_(trImg, nFrames_led)
if nargin<2, nFrames_led = []; end
if ~isempty(nFrames_led)
    nFrames_led = min(size(trImg,3), nFrames_led);
    trImg = trImg(:,:,1:nFrames_led);
end
img_pp = (max(trImg,[],3) - min(trImg,[],3));
[~,imax_pp] = max(img_pp(:));
[yLed, xLed] = ind2sub(size(img_pp), imax_pp);
xyLed = [xLed, yLed];
end %func


%--------------------------------------------------------------------------
% 7/30/2018 JJJ: Moved from GUI.m
function vcFile_Track = trial_save_(handles)

handles.ESAC = calcESAC(handles);

[handles.vcVer, handles.vcVer_date] = version_();
S_cfg = vistrack('load-cfg');
S_save = struct_copy_(handles, S_cfg.csFields);

if isfield(handles, 'vcFile_Track')
    vcFile_Track = handles.vcFile_Track;
elseif exist_file_(handles.vidFname)
    vcFile_Track = subsFileExt_(handles.vidFname, '_Track.mat');
else
    vcFile_Track = get(handles.editResultFile, 'String');
end

h = msgbox('Saving... (this will close automatically)');    
try
    struct_save_(S_save, vcFile_Track, 0);
    if isfield(handles, 'editResultFile')
        set(handles.editResultFile, 'String', vcFile_Track);
        msgbox_(sprintf('Output saved to %s', fullpath_(vcFile_Track)));   
    else
        fprintf('Output saved to %s\n', fullpath_(vcFile_Track)); % batch mode
    end
catch
    fprintf(2, 'Save file failed: %s\n', vcFile_Track);
end
close_(h);
end %func


%--------------------------------------------------------------------------
% 7/30/2018 JJJ: Moved from GUI.m
function S_save = struct_copy_(handles, csField)
for i=1:numel(csField)
    try
        S_save.(csField{i}) = handles.(csField{i});
    catch
        S_save.(csField{i}) = []; % not copied
    end
end
end %func


%--------------------------------------------------------------------------
% 7/24/2018: Copied from jrc3.m
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end
end %func


%--------------------------------------------------------------------------
% 7/30/18 JJJ: Copied from jrc3.m
function struct_save_(S, vcFile, fVerbose)
nRetry = 3;
if nargin<3, fVerbose = 0; end
if fVerbose
    fprintf('Saving a struct to %s...\n', vcFile); t1=tic;
end
version_year = version('-release');
version_year = str2double(version_year(1:end-1));
if version_year >= 2017
    for iRetry=1:nRetry
        try
            save(vcFile, '-struct', 'S', '-v7.3', '-nocompression'); %faster    
            break;
        catch
            pause(.5);
        end
        fprintf(2, 'Saving failed: %s\n', vcFile);
    end
else    
    for iRetry=1:nRetry
        try
            save(vcFile, '-struct', 'S', '-v7.3');   
            break;
        catch
            pause(.5);
        end
        fprintf(2, 'Saving failed: %s\n', vcFile);
    end    
end
if fVerbose
    fprintf('\ttook %0.1fs.\n', toc(t1));
end
end %func


%--------------------------------------------------------------------------
function [S_sync, mov] = calc_sync_(handles, mov)
% handles.{ADCTS, vidFname, vidobj}

if nargin<2, mov = []; end
if isempty(mov), mov = handles2mov_(handles); end

% Find LED timing
vtLed_adc = getSync_adc_(handles);
[vtLed_adc, viLed_adc_removed] = remove_pulses_(vtLed_adc);

xyLed = findLed_mov_(mov, 300);
vrLed = mov2led_(mov, xyLed);
vrLed = vrLed - medfilt1(vrLed,5);
thresh_led = max(vrLed) * .2;
viLed_cam = find(diff(vrLed > thresh_led)>0) + 1;
[viLed_cam, viiLed_filled] = fill_pulses_(viLed_cam);

if numel(viLed_cam) > numel(vtLed_adc), viLed_cam(1) = []; end % remove the first
S_sync = struct('vrT_adc', vtLed_adc, 'viT_cam', viLed_cam);
end %func


%--------------------------------------------------------------------------
function mov = handles2mov_(handles, P)
if nargin<2, P = []; end
if isempty(P), P = load_settings_(handles); end

vcFile_vid = handles.vidFname;
vcVidExt = get_set_(P, 'vcVidExt', '.wmv');
if ~exist_file_(vcFile_vid)
    vcFile_vid = strrep(get_(handles, 'vcFile_Track'), '_Track.mat', vcVidExt);
end
h = msgbox_('Loading video... (this closes automatically)');
[mov, dimm_vid] = loadvid_(vcFile_vid, get_set_(P, 'nSkip_vid', 4));
close_(h);
end %func


%--------------------------------------------------------------------------
function [handles, hFig] = trial_fixsync_(handles, fPlot)
% fPlot: 0 (no-plot, save), 1 (plot, save), 2 (plot, no save)
persistent mov

if nargin==0, mov = []; return; end % clear cache
if nargin<2, fPlot = 1; end

P = load_settings_(handles);
if isempty(mov)
    mov = handles2mov_(handles, P);
else
    fprintf('Using cached video.\n');
end

S_sync = calc_sync_(handles, mov);
% [vrT_adc, viT_cam] = struct_get_(S_sync, 'vrT_adc', 'viT_cam');
[tlim_adc, flim_cam] = sync_limit_(S_sync.vrT_adc, S_sync.viT_cam);

% save if not plotting
if fPlot == 0
    TC = cam2adc_sync_(S_sync, handles.FLIM(1):handles.FLIM(2));
    FPS = diff(handles.FLIM([1,end])) / diff(handles.TC([1,end]));    
    save(handles.vcFile_Track, 'TC', 'FPS', 'S_sync', '-append');
    return;
end

[vtText, csText] = getText_adc_(handles, P);
xoff_ = 50;
csPopup = {'First frame', csText{:}, 'Last frame'};
hFig = figure_new_('FigSync', [handles.vidFname, ' press "h" for help'], [0,0,.5,1]);
hFig.KeyPressFcn = @keypress_FigSync_;
hPopup = uicontrol('Style', 'popup', 'String', csPopup, ...
    'Position', [xoff_ 0 200 50], 'Callback', @popup_sync_);
hPopup.KeyPressFcn = @(h,e)keypress_FigSync_(hFig,e); 

% Create axes
iFrame = 1;
hAxes1 = axes(hFig, 'Units', 'pixels', 'Position', [xoff_,60,800,600]);
hImage = imshow(mov(:,:,iFrame), 'Parent', hAxes1);
hold(hAxes1, 'on');
hTitle = title_(hAxes1, sprintf('Frame %d', iFrame));

% Create Line plot
tFrame_adc = cam2adc_sync_(S_sync, iFrame);
hAxes2 = axes(hFig, 'Units', 'pixels', 'Position', [xoff_,750,800,50]);
hold(hAxes2, 'on');
plot_vline_(hAxes2, S_sync.vrT_adc, [0,1], 'k'); 
plot_vline_(hAxes2, vtText, [0,1], 'm'); 
hLine_cam = plot_vline_(hAxes2, cam2adc_sync_(S_sync, S_sync.viT_cam), [0,1], 'r--'); 
set(hAxes2, 'XTick', S_sync.vrT_adc);
xlabel('ADC Time (s)');
set(hAxes2, 'XLim', tlim_adc);
hCursor_adc = plot(hAxes2, tFrame_adc, .5, 'ro');

% Create Line plot
hAxes3 = axes(hFig, 'Units', 'pixels', 'Position', [xoff_,850,800,50]);
hold(hAxes3, 'on');
hPlot3 = plot_vline_(hAxes3, S_sync.viT_cam, [0,1], 'r--');
xlabel('Camera Frame #');
set(hAxes3, 'XLim', flim_cam);
set(hAxes3, 'XTick', S_sync.viT_cam);
hCursor_cam = plot(hAxes3, iFrame, .5, 'ro');

% Show EOD plot
hAxes4 = axes(hFig, 'Units', 'pixels', 'Position', [xoff_,950,800,100]);
hold(hAxes4, 'on');
xlabel('ADC Time (s)');
ylabel('EOD Rate (Hz)');
[vtEodr, vrEodr] = getEodr_adc_(handles, P);
hPlot_eod = plot(hAxes4, vtEodr, vrEodr, 'k');
hCursor_eod = plot(hAxes4, tFrame_adc, median(vrEodr), 'ro');
set(hAxes4, 'XLim', tlim_adc, 'YLim', median(vrEodr) * [1/2, 2]);

hFig.UserData = makeStruct_(iFrame, mov, hImage, vtText, csText, hTitle, ...
    S_sync, hCursor_adc, hCursor_cam, hPopup, hPlot_eod, hCursor_eod);
set0_(S_sync);

if fPlot == 2, return; end

% close the figure after done
msgbox_('Close the figure when finished.');
uiwait(hFig);

S_sync = get0_('S_sync');
TC = cam2adc_sync_(S_sync, handles.FLIM(1):handles.FLIM(2));    
if questdlg_(sprintf('Save time sync? (mean error: %0.3fs)', std(TC-handles.TC)))
    handles.TC = TC;
    handles.FPS = diff(handles.FLIM([1,end])) / diff(handles.TC([1,end]));
    trial_save_(handles);
end
end %func


%--------------------------------------------------------------------------
function vrT1_adc = cam2adc_sync_(S_sync, viT1_cam)
if isempty(S_sync), S_sync = get0_('S_sync'); end
[vrT_adc, viT_cam] = struct_get_(S_sync, 'vrT_adc', 'viT_cam');
nBlinks = min(numel(viT_cam), numel(vrT_adc));
[viT_cam, vrT_adc] = deal(viT_cam(end-nBlinks+1:end), vrT_adc(1:nBlinks));
vrT1_adc = interp1(viT_cam, vrT_adc, viT1_cam, 'linear', 'extrap');
end %func


%--------------------------------------------------------------------------
function [tlim_adc, flim_cam] = sync_limit_(vtLed_adc, viLed_cam)
nBlinks = min(numel(vtLed_adc), numel(viLed_cam));
tlim_adc = vtLed_adc([1, nBlinks]);
flim_cam = viLed_cam([end-nBlinks+1, end]);
end %func


%--------------------------------------------------------------------------
function viT1_cam = adc2cam_sync_(S_sync, vrT1_adc)
if isempty(S_sync), S_sync = get0_('S_sync'); end
[vrT_adc, viT_cam] = struct_get_(S_sync, 'vrT_adc', 'viT_cam');
nBlinks = min(numel(viT_cam), numel(vrT_adc));
[viT_cam, vrT_adc] = deal(viT_cam(1:nBlinks), vrT_adc(1:nBlinks));
viT1_cam = round(interp1(vrT_adc, viT_cam, vrT1_adc, 'linear', 'extrap'));
end %func


%--------------------------------------------------------------------------
function vc = popup_sync_(h,e)
hFig = h.Parent;
S_fig = hFig.UserData;
vcLabel = h.String{h.Value};
[iFrame_prev, mov, S_sync, hTitle, hImage] = ...
    struct_get_(S_fig, 'iFrame', 'mov', 'S_sync', 'hTitle', 'hImage');
nFrames = size(mov,3);
switch lower(vcLabel)
    case 'first frame'
        iFrame = 1;
    case 'last frame'
        iFrame = nFrames; 
    otherwise
        t_adc = S_fig.vtText(h.Value-1);
        iFrame = setlim_(adc2cam_sync_(S_sync, t_adc), [1, nFrames]);
end
if iFrame_prev==iFrame, return ;end

refresh_FigSync_(hFig, iFrame);
end %func


%--------------------------------------------------------------------------
function refresh_FigSync_(hFig, iFrame)
S_fig = hFig.UserData;
[iFrame_prev, mov, S_sync, hTitle, hImage] = ...
    struct_get_(S_fig, 'iFrame', 'mov', 'S_sync', 'hTitle', 'hImage');
hImage.CData = mov(:,:,iFrame);
set(S_fig.hCursor_cam, 'XData', iFrame);
tFrame_adc = cam2adc_sync_(S_sync, iFrame);
set(S_fig.hCursor_adc, 'XData', tFrame_adc);
set(S_fig.hCursor_eod, 'XData', tFrame_adc);

% Update title
viText_cam = adc2cam_sync_(S_sync, S_fig.vtText);
viMatch = find(viText_cam==iFrame);
if isempty(viMatch)
    hTitle.String = sprintf('Frame %d', iFrame);
else
    iMatch = viMatch(1);
    hTitle.String = sprintf('Frame %d: %s', iFrame, S_fig.csText{iMatch});
    S_fig.hPopup.Value = iMatch+1;
end

% update current frame
S_fig.iFrame = iFrame;
hFig.UserData = S_fig;
set0_(S_fig); % push to global
end %func


%--------------------------------------------------------------------------
function vr = setlim_(vr, lim_)
% Set low and high limits
vr = min(max(vr, lim_(1)), lim_(2));
end %func


%--------------------------------------------------------------------------
function hTitle = title_(hAx, vc)
% title_(vc)
% title_(hAx, vc)

if nargin==1, vc=hAx; hAx=[]; end
% Set figure title

if isempty(hAx), hAx = gca; end
hTitle = get_(hAx, 'Title');
if isempty(hTitle)
    hTitle = title(hAx, vc, 'Interpreter', 'none', 'FontWeight', 'normal');
else
    set_(hTitle, 'String', vc, 'Interpreter', 'none', 'FontWeight', 'normal');
end
end %func


%--------------------------------------------------------------------------
function vc = set_(vc, varargin)
% Set handle to certain values
% set_(S, name1, val1, name2, val2)

if isempty(vc), return; end
if isstruct(vc)
    for i=1:2:numel(varargin)        
        vc.(varargin{i}) = varargin{i+1};
    end
    return;
end
if iscell(vc)
    for i=1:numel(vc)
        try
            set(vc{i}, varargin{:});
        catch
        end
    end
elseif numel(vc)>1
    for i=1:numel(vc)
        try
            set(vc(i), varargin{:});
        catch
        end
    end
else
    try
        set(vc, varargin{:});
    catch
    end 
end
end %func


%--------------------------------------------------------------------------
function clear_cache_()
trial_fixsync_();
end %func


%--------------------------------------------------------------------------
function [mov, vcFile_bin] = loadvid_(vcFile_vid, nSkip, fSave_bin)
% using the 2018a VideoReader
% Extracts red channel only
t1=tic;
if nargin<2, nSkip = []; end
if nargin<3, fSave_bin = []; end

if isempty(nSkip), nSkip = 1; end
if isempty(fSave_bin), fSave_bin = 1; end

fFast = 0; %subsampling instead of averaging the pixels

try
    vidobj = VideoReader(vcFile_vid);
catch
    [dimm, mov] = deal([]);
    return;
end

% vidInfo = mmfileinfo(vcFile_vid);
% vidInfo.Duration;

fprintf('Loading video: %s\n', vcFile_vid);

vidHeight = floor(vidobj.Height / nSkip);
vidWidth = floor(vidobj.Width / nSkip);

% check cache
vcFile_bin = sprintf('%s_mov%dx%d.bin', vcFile_vid, vidHeight, vidWidth);
mov = loadvid_bin_(vcFile_bin);
if ~isempty(mov)
    dimm = size(mov); 
    fprintf('\tLoaded from %s (%d frames), took %0.1fs\n', vcFile_bin, size(mov,3), toc(t1));
    return; 
end

nFrames_est = round(vidobj.Duration * vidobj.FrameRate);
mov = zeros(vidHeight, vidWidth, nFrames_est, 'uint8'); 
fTrim = (vidHeight * nSkip) < vidobj.Height || (vidWidth * nSkip) < vidobj.Width;
iFrame = 0;
while hasFrame(vidobj)
    iFrame = iFrame + 1;    
    img_ = readFrame(vidobj);
    if fFast
        mov(:,:,iFrame) = img_(1:nSkip:vidHeight*nSkip, 1:nSkip:vidWidth*nSkip, 1);
        continue;
    end
    img_ = img_(:,:,1); % red extraction
    if nSkip>1
        if fTrim
            img_ = img_(1:(vidHeight*nSkip), 1:(vidWidth*nSkip));
        end
        img_ = sum(uint16(reshape(img_, nSkip, [])));
        img_ = sum(permute(reshape(img_, vidHeight, nSkip, vidWidth), [2,1,3]));
        img_ = reshape(uint8(img_/(nSkip^2)), vidHeight, vidWidth);
    end
    mov(:,:,iFrame) = img_;
end
nFrames = iFrame;
dimm = [vidHeight, vidWidth, nFrames];
if nFrames < nFrames_est
    mov = mov(:,:,1:nFrames); %trim
end

% bulk save
if fSave_bin
    try 
        fid_w = fopen(vcFile_bin, 'w');
        fwrite(fid_w, mov, class(mov));
        fclose(fid_w); 
        fprintf('\twrote to %s (%d frames), took %0.1fs\n', vcFile_bin, size(mov,3), toc(t1));
    catch
        ;
    end    
else
    fprintf('\tLoaded %d frames, took %0.1fs\n', size(mov,3), toc(t1));
end
end %func


%--------------------------------------------------------------------------
function mov = loadvid_bin_(vcFile_bin)
% vcFile_bin: string format: vidfile_mov%dx%d.bin (wxh)

mov=[]; 
if ~exist_file_(vcFile_bin), return; end
vcFormat = regexpi(vcFile_bin, '_mov(\d+)[x](\d+)[.]bin$', 'match');
if isempty(vcFormat), return; end % invalid format

try
    vcFormat = strrep(strrep(vcFormat{1}, '_mov', ''), '.bin', '');
    dimm = sscanf(vcFormat, '%dx%d');
    [height, width] = deal(dimm(1), dimm(2));
    nBytes_file = filesize_(vcFile_bin);
    dimm(3) = floor(nBytes_file/height/width);
    
    fid = fopen(vcFile_bin, 'r');
    mov = fread_(fid, dimm, 'uint8');
    fclose(fid);
catch
    return;
end
end %func


%--------------------------------------------------------------------------
function mnWav1 = fread_(fid_bin, dimm_wav, vcDataType)
% Get around fread bug (matlab) where built-in fread resize doesn't work
dimm_wav = dimm_wav(:)';
try
    if isempty(dimm_wav)
        mnWav1 = fread(fid_bin, inf, ['*', vcDataType]);
    else
        if numel(dimm_wav)==1, dimm_wav = [dimm_wav, 1]; end
        mnWav1 = fread(fid_bin, prod(dimm_wav), ['*', vcDataType]);
        if numel(mnWav1) == prod(dimm_wav)
            mnWav1 = reshape(mnWav1, dimm_wav);
        else
            dimm2 = floor(numel(mnWav1) / dimm_wav(1));
            if dimm2 >= 1
                mnWav1 = reshape(mnWav1, dimm_wav(1), dimm2);
            else
                mnWav1 = [];
            end
        end
    end
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
% Return [] if multiple files are found
function nBytes = filesize_(vcFile)
S_dir = dir(vcFile);
if numel(S_dir) ~= 1
    nBytes = []; 
else
    nBytes = S_dir(1).bytes;
end
end %func



%--------------------------------------------------------------------------
function [S_trialset, trFps] = trialset_fixfps_(vcFile_trialset)
% It loads the files
% iData: 1, ang: -0.946 deg, pixpercm: 7.252, x0: 793.2, y0: 599.2
% run S141106_LearningCurve_Control.m first cell
fFix_sync = 1;

S_trialset = load_trialset_(vcFile_trialset);
% [pixpercm, angXaxis] = struct_get_(S_trialset.P, 'pixpercm', 'angXaxis');
[tiImg, vcType_uniq, vcAnimal_uniq, viImg, csFiles_Track] = ...
    struct_get_(S_trialset, 'tiImg', 'vcType_uniq', 'vcAnimal_uniq', 'viImg', 'csFiles_Track');

hMsg = msgbox('Analyzing... (This closes automatically)');
t1=tic;
trFps = nan(size(tiImg));
for iTrial = 1:numel(viImg)    
    try
        clear_cache_();
        S_ = load(csFiles_Track{iTrial}, 'TC', 'XC', 'YC', 'xy0', 'vidFname', 'FPS', 'img0', 'ADCTS', 'FLIM');         
        S_.vcFile_Track = csFiles_Track{iTrial};
        if fFix_sync, S_ = trial_fixsync_(S_, 0); end
        iImg_ = viImg(iTrial);        
        trFps(iImg_) = get_set_(S_, 'FPS', nan);     
        fprintf('\n');
    catch
        disp(csFiles_Track{iTrial});
    end
end %for
fprintf('\n\ttook %0.1fs\n', toc(t1));
close_(hMsg);

if nargout==0
    hFig = plot_trialset_img_(S_trialset, trFps); 
    set(hFig, 'Name', sprintf('FPS: %s', vcFile_trialset));
end
end %func


%--------------------------------------------------------------------------
% 8/9/2018 JJJ: copied from irc.m
function varargout = get0_(varargin)
% returns get(0, 'UserData') to the workspace
% [S0, P] = get0_();
S0 = get(0, 'UserData'); 
if nargin==0
    if nargout==0
        assignWorkspace_(S0); 
    else
        varargout{1} = S0; 
    end
else
    for iArg=1:nargin
        try                
            eval(sprintf('%s = S0.%s;', varargin{iArg}, varargin{iArg}));
            varargout{iArg} = S0.(varargin{iArg});
        catch
            varargout{iArg} = [];
        end
    end
end
end %func


%--------------------------------------------------------------------------
% 8/9/2018 JJJ: copied from irc.m
function S0 = set0_(varargin)

S0 = get(0, 'UserData'); 
for i=1:nargin
    try
        S0.(inputname(i)) = varargin{i};
    catch
        disperr_();
    end
end
set(0, 'UserData', S0);
end %func


%--------------------------------------------------------------------------
function mov = loadvid_preview_(vcFile_vid, viFrames)
if nargin<2, viFrames = []; end

if ~ischar(vcFile_vid)
    vcFile_vid = fullfile(vcFile_vid.Path, vcFile_vid.Name);
end
P = load_cfg_();
mov = loadvid_(vcFile_vid, get_set_(P, 'nSkip_vid', 4));
if ~isempty(viFrames), mov = mov(:,:,viFrames); end
end %func


%--------------------------------------------------------------------------
function handles = trial_sync_(handles)
[handles.S_sync, mov] = calc_sync_(handles);
[vrT_adc, viT_cam] = struct_get_(handles.S_sync, 'vrT_adc', 'viT_cam');
handles.TLIM0 = vrT_adc([1, end]);
handles.FLIM0 = viT_cam([1, end]);
handles.FPS = diff(handles.FLIM0) / diff(handles.TLIM0);

% plot sync
[~, hFig] = trial_fixsync_(handles, 2);
msgbox({'Close the figure after checking the sync.', 'Press PageUp/PageDown/Left/Right to navigate'});
uiwait(hFig);

vcAns = questdlg('Synchronized correctly?');
if strcmpi(vcAns, 'Yes')
    set(handles.btnBackground, 'Enable', 'on');
else
    set(handles.btnBackground, 'Enable', 'off');
end
end %func


%--------------------------------------------------------------------------
function trialset_import_track_(vcFile_trialset)

% Find destination
S_trialset = load_trialset_(vcFile_trialset);
if isempty(S_trialset), errordlg('No trials exist', vcFile_trialset); return; end
vcVidExt = get_set_(S_trialset.P, 'vcVidExt');
[csFiles_vid, csDir_vid] = find_files_(S_trialset.vcDir, ['*', vcVidExt]);

% ask from where
vcDir_copyfrom = fileparts(S_trialset.vcDir);
vcDir_copyfrom = uigetdir(vcDir_copyfrom, 'Select a folder to copy from');
if ~ischar(vcDir_copyfrom), return; end
csFiles_Track = find_files_(vcDir_copyfrom, '*_Track.mat');
if isempty(csFiles_Track), return; end

fprintf('Copying %d files\n', numel(csFiles_Track));
nCopied = 0;
for iFile_Track = 1:numel(csFiles_Track)
    try
        vcFile_from_ = csFiles_Track{iFile_Track};

        [~,vcFile_to_,~] = fileparts(vcFile_from_);
        vcFile_to_ = cellstr_find_(csFiles_vid, strrep(vcFile_to_, '_Track', vcVidExt));
        vcFile_to_ = strrep(vcFile_to_, vcVidExt, '_Track.mat');
        
        copyfile(vcFile_from_, vcFile_to_, 'f');        
        fprintf('\tCopying %s to %s\n', vcFile_from_, vcFile_to_);
        nCopied = nCopied + 1;
    catch
        fprintf(2, '\tCopy error: %s to %s\n', vcFile_from_, vcFile_to_);
    end
end %for
fprintf('\t%d/%d copied\n', nCopied, numel(csFiles_Track));
end %func


%--------------------------------------------------------------------------
function vc_match = cellstr_find_(csFrom, vcFind)
cs = cellfun(@(vcFrom)regexpi(vcFrom, vcFind, 'match'), csFrom, 'UniformOutput', 0);
iFind = find(~cellfun(@isempty, cs));
if isempty(iFind)
    vc_match = []; 
else
    vc_match = csFrom{iFind(1)};
end
end %func


%--------------------------------------------------------------------------
function trialset_googlesheet_(vcFile_trialset)
S_trialset = load_trialset_(vcFile_trialset);
vcLink_googlesheet = get_(S_trialset, 'vcLink_googlesheet');
if isempty(vcLink_googlesheet)
    fprintf('"vcLink_googlesheet" is not set in %d\n', vcArg1);
else
    web_(vcLink_googlesheet);
end
end %func


%--------------------------------------------------------------------------
function prefix = getSpike2Prefix_(S)
    prefix = fields(S); 
    prefix = prefix{1}; 
    k = strfind(prefix, '_Ch');
    k=k(end);
    prefix = prefix(1:k-1);    
end