function plotAnimals_EODA(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, strVar, fun1)
% plot correlatoin coefficient

csPhase = {'E', 'L', 'P'};
csAnimal = {'A', 'B', 'C', 'D', 'All'};
csZone = {'AZ', 'LM', 'NF', 'F'};
cvLM = cell(4,3);
mrLM = zeros(4,3);

figure;
% suptitle([strVar ', ' func2str(fun1)]);

%-------------------
% Plot per animal stats
for iAnimal = 1:numel(csAnimal)
    subplot(1,5,iAnimal);   
    for iZone = 1:numel(csZone);        
        for iPhase = 1:numel(csPhase)
            eval(sprintf('vsTrialPool = vsTrialPool_%s;', csPhase{iPhase}));
    %         S = poolTrials_location(vsTrialPool, iAnimal);
            if iAnimal <= 4
                S = poolTrials_IPI(vsTrialPool, iAnimal);            
            else
                S = poolTrials_IPI(vsTrialPool, []); 
            end
            [vlZone, strZone] = getZone(S, iZone);
            eval(sprintf('vrZ = %s;', strVar));
            mrLM(iZone, iPhase) = fun1(vrZ(vlZone));
        end        
    end
    h = bar(mrLM);
    set(gca, 'XTickLabel', csZone);    
    set(h(1), 'FaceColor', 'r');
    set(h(2), 'FaceColor', 'b');
    set(h(3), 'FaceColor', 'g');    
    title(csAnimal{iAnimal});
end %for

end %func


function [vlZone, strZone] = getZone(S, iZone)
switch (iZone)
    case 1 %all
        vlZone = S.vlZone;
        strZone = 'AZ';
    case 2 %LM
        vlZone = S.vrD1 <= 3 | S.vrD2 <= 3 | S.vrD3 <= 3 | S.vrD4 <= 3; %within landmark detection zone
        strZone = 'LM<3';
    case 3 %Fc<15
        vlZone = S.vrDf < 14 & S.vrDf >= 3;
        strZone = 'Fc4~15';
    case 4 %F<3
        vlZone = S.vrDf < 3;
        strZone = 'F<3';
end
end %func