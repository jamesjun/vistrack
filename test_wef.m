S_out = wef('test', 'mfile2struct_', {'settings_wef.m'}, 1)

vcFile_mat = 'G:\Malerlab\Tracked spatial learning\2013 E-complete\E13D3_Track.mat';
S_out = wef('test', 'checkFile_', {vcFile_mat, 'D'}, 2);
S_out = wef('test', 'checkFile_', {vcFile_mat, 'A'}, 2);
vcFile_mat = 'G:\Malerlab\Tracked spatial learning\2013 E-complete\E14B3p_Track.mat';
S_out = wef('test', 'checkFile_', {vcFile_mat, 'B'}, 2);
S_out = wef('test', 'checkFile_', {vcFile_mat, 'A'}, 2);

S_out = wef('test', 'vc2cs_', {'A'}, 1);
S_out = wef('test', 'vc2cs_', {'AB'}, 1);
S_out = wef('test', 'vc2cs_', {'ABC'}, 1);
S_out = wef('test', 'vc2cs_', {''}, 1);
