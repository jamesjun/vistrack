function csZ = plotPhase(csTrials, iZone, vcVar, iAnimal)
   
if nargin < 4
    iAnimal = [];
end

%csTrials = {vsTrialPool_E, vsTrialPool_L, vsTrialPool_P}; 
csXlabel = {'Early', 'Late', 'Probe'};        
 
csZ = cell(csTrials); 
for iPhase = 1:numel(csTrials)
    vsTrials = csTrials{iPhase};
    if isempty(iAnimal)
        viTrials = 1:numel(vsTrials);
    else
        viTrials = find([vsTrials.iAnimal] == iAnimal);
    end
    vrZ = nan(size(viTrials));
    for iTrial1 = 1:numel(viTrials)
        iTrial = viTrials(iTrial1);
        RS = poolTrials_RS(vsTrials(iTrial), iAnimal, [], iZone);   
        eval(sprintf('z = RS.%s;', vcVar));
        z = z(RS.vlZ0);
        if ~isempty(z)
            vrZ(iTrial1) = nanmean(z);
        end
    end
    csZ{iPhase} = vrZ;
end

if nargout==0
    boxplot_cell(csZ, 'mean-sem'); %median +- iqr
    %[pval, k, K] = circ_kuipertest(csAngErr{1}/90*2*pi, csAngErr{2}/90*2*pi)
    set(gca, 'XTickLabel', csXlabel);
end