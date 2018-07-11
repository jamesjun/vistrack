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
    otherwise, help_(); return;
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
csFiles_upload = {'*.m', 'GUI.fig', 'change_log.txt', 'readme.txt'};
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
vcVer = 'v0.1.2';
vcDate = '7/11/2018';
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
    '[Deployment]';
    '  vistrack update';
    '    Update from Github'; 
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


