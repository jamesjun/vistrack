function plotAnimals_var(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, strVar, fun1)
% plot correlatoin coefficient

vsPhase = {'E', 'L', 'P'};
cvLM = cell(4,3);
mrLM = zeros(4,3);

figure;
suptitle([strVar ', ' func2str(fun1)]);

%-------------------
% Plot per animal stats
for iZone = 1:4
    subplot(3,2,iZone);
    for iAnimal = 1:4
        for iPhase = 1:3
            eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    %         S = poolTrials_location(vsTrialPool, iAnimal);
            S = poolTrials_IPI(vsTrialPool, iAnimal);            
            [vlZone, strZone] = getZone(S, iZone);
            eval(sprintf('vrZ = S.%s;', strVar));
            mrLM(iAnimal, iPhase) = fun1(vrZ(vlZone));
        end    
    end
    h = bar(mrLM);
    set(h(1), 'FaceColor', 'r');
    set(h(2), 'FaceColor', 'b');
    set(h(3), 'FaceColor', 'g');
    set(gca, 'XTickLabel', {'A', 'B', 'C', 'D'});    
    title(strZone);
end %for

%-------------------
% Plot per loc stats
subplot(3,2,5);
mrLM = zeros(4,3);
for iZone = 1:4
    for iPhase = 1:3
        eval(sprintf('vsTrialPool = vsTrialPool_%s;', vsPhase{iPhase}));
    %         S = poolTrials_location(vsTrialPool, iAnimal);
        S = poolTrials_IPI(vsTrialPool, iAnimal);            
        [vlZone, strZone] = getZone(S, iZone);
        eval(sprintf('vrZ = S.%s;', strVar));
        mrLM(iZone, iPhase) = fun1(vrZ(vlZone));
    end      
end %for
h = bar(mrLM);
title('All animals');
set(h(1), 'FaceColor', 'r');
set(h(2), 'FaceColor', 'b');
set(h(3), 'FaceColor', 'g');  
set(gca, 'XTickLabel', {'AZ', 'LM<3', 'Fc<15', 'F<3'});    

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