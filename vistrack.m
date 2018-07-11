function varargout = vistrack(varargin)
vcCmd = 'help';
if nargin==1, vcCmd = varargin{1}; end

switch vcCmd
    case 'commit',  commit_(); return;
    case 'help',
end %switch

end %func


%--------------------------------------------------------------------------
function commit_(vcDir_target)
if nargin<1, vcDir_target = 'D:\Dropbox\Git\vistrack\'; end

delete_empty_files_();
delete([vcDir_target, '*']);
copyfile('*.m', vcDir_target, 'f');
copyfile('GUI.fig', vcDir_target, 'f');
copyfile('change_log.txt', vcDir_target, 'f');
copyfile('readme.txt', vcDir_target, 'f');
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