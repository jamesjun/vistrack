function csZ = plotPhaseIpi(csTrials, iZone, vcVar, iAnimal)
   
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
        IPI = poolTrials_IPI(vsTrials(iTrial), iAnimal, [], iZone);   
        eval(sprintf('z = %s;', vcVar));
        vrZ(iTrial1) = z;
    end
    csZ{iPhase} = vrZ;
end

if nargout==0
    boxplot_cell(csZ, 'mean-sem'); %median +- iqr
    %[pval, k, K] = circ_kuipertest(csAngErr{1}/90*2*pi, csAngErr{2}/90*2*pi)
    set(gca, 'XTickLabel', csXlabel);
    xlim([.5 3.5]);
end