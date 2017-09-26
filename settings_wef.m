%--------------------------------------------------------------------------
% common setting
angXaxis = -1.1590;
pixpercm = 1053.28/(sqrt(2)*100);
csNames_set = {'rand', 'randwide', 'none', 'cue', 'stable', 'shuffle'};
nTrialsPerSession = 4;
quantLim_dist = [1/8, 7/8];
maxDur = 360; % in sec
csAnimals = {'A','B','C','D'};
rectCrop = [493 1083 312 902];

%--------------------------------------------------------------------------
% Set-specific info

vcDir_rand = 'G:\Malerlab\Tracked spatial learning\2013b_RandomNarrow_complete\';
csFiles_rand = 'G:\Malerlab\2013b\*_R*\';
% vcDataset_rand = 'D141026_RandGroup.mat';
vcDataset_rand = 'Dataset_rand.mat';
csAnimals_rand = csAnimals;

vcDir_randwide = 'G:\Malerlab\Tracked spatial learning\2013b_WideRandom_complete\';
csFiles_randwide = 'G:\Malerlab\2013b\*_W*\';
% vcDataset_randwide = 'D141026_RandWideGroup.mat';
vcDataset_randwide = 'Dataset_randwide.mat';
csAnimals_randwide = csAnimals;

vcDir_none = 'G:\Malerlab\Tracked spatial learning\2013b_None_complete\';
csFiles_none = 'G:\Malerlab\2013b\*_N*\';
% vcDataset_none = 'D141026_NoneGroup';
vcDataset_none = 'Dataset_none.mat';
csAnimals_none = csAnimals;

vcDir_stable = 'G:\Malerlab\Tracked spatial learning\2013 E-complete\';
csFiles_stable = 'G:\Malerlab\2013b\*_E*\';
% vcDataset_stable = 'D140324_LandmarkGroup.mat';
vcDataset_stable = 'Dataset_stable.mat';
csAnimals_stable = csAnimals;
cviProbe_stable = {[13:16], [17:20], [21:24], [25:26]};
csProbe_stable = {'probe', 'contraction', 'expansion', 'translation'};

vcDir_relearn = 'G:\Malerlab\Tracked spatial learning\2013 E-complete\';
csFiles_relearn = 'G:\Malerlab\2013b\*_E*\';
% vcDataset_stable = 'D140324_LandmarkGroup.mat';
vcDataset_relearn = 'Dataset_relearn.mat';
csAnimals_relearn = csAnimals;

vcDir_cue = 'G:\Malerlab\Tracked spatial learning\2013b_Cue_complete\';
csFiles_cue = 'G:\Malerlab\2013b\*_C*\';
vcDataset_cue = 'Dataset_cue.mat';
csAnimals_cue = csAnimals;

vcDir_shuffle = 'G:\Malerlab\Tracked spatial learning\2014b_ShuffleCue\';
csFiles_shuffle = 'G:\Malerlab\2014b\*_S*\';
vcDataset_shuffle = 'Dataset_shuffle.mat';
csAnimals_shuffle = csAnimals;
% fish A,B: stable cue (large circle, large square), C,D: random cue


