S_trialset = load('C:\Users\James Jun\Downloads\cue_trialset.mat');
vS_trial=cell2mat(S_trialset.cS_trial);
cmrPos_shape={vS_trial.mrPos_shape};
[~,csDataID] = cellfun(@fileparts, {vS_trial.vidFname}, 'UniformOutput', 0);
viAnimal = cellfun(@(x)x(4)-'A'+1, csDataID);

vlFilled = cellfun(@(x)all(any(isnan(x),2)), cmrPos_shape)

mrPos_all = cell2mat(cellfun(@(x)x(:), cmrPos_shape, 'UniformOutput', 0));
% vrUniq = unique(mrPos_all(:));
% vrUniq = vrUniq(~isnan(vrUniq)); % nan is represented separately
% miPos_all = ismember(mrPos_all, vrUniq);
mrPos_all(isnan(mrPos_all)) = nan;
[C,ia,ic] = unique(mrPos_all', 'rows') % 21 types found?

vrPos_unique = unique(mrPos_all(:))
mrPos_(isnan(mrPos_)) = 2^15;
mrPos_ = int16(mrPos_);




%% find empty 
