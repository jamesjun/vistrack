function [AX, AX1, cvZ] = plotAnimals(vsTrialPool_E, vsTrialPool_L, vsTrialPool_P, strVar, fun1, strY)
% cvZ: cell of animal, zone, phase
% plot correlatoin coefficient

csPhase = {'E', 'L', 'P'};
csPhaseColor = {'r', 'b', 'g'};
csDiffColor = {'m', 'c'};
csAnimal = {'All', 'A', 'B', 'C', 'D'};
csZone = {'AZ', 'LM', 'NF', 'F'};
mrMean = zeros(4,3);
mrSem = zeros(4,3);

vlDep = @(S)(S.vrD1 <= 3 & differentiate3(S.vrD1) > 0) |...
               (S.vrD2 <= 3 & differentiate3(S.vrD2) > 0) |...
               (S.vrD3 <= 3 & differentiate3(S.vrD3) > 0) |...
               (S.vrD4 <= 3 & differentiate3(S.vrD4) > 0);
         
vlApp = @(S)(S.vrD1 <= 3 & differentiate3(S.vrD1) < 0) |...
               (S.vrD2 <= 3 & differentiate3(S.vrD2) < 0) |...
               (S.vrD3 <= 3 & differentiate3(S.vrD3) < 0) |...
               (S.vrD4 <= 3 & differentiate3(S.vrD4) < 0);
           
if nargin < 5
    fun1 = [];
end
if isempty(fun1)
    fun1 = @(x)nanmean(x);
%     fun2 = @(x)nanstd(x);
    fun2 = @(x)nanstd(x) / sqrt(numel(x)); %sem
else
    fun2 = @(x)0;
end

if nargin < 6
    strY = strVar;
end

figure;
% suptitle([strVar ', ' func2str(fun1)]);

%-------------------
% Plot per animal stats
AX = []; %individual bars
AX1 = []; %change bars
cvZ = cell(numel(csAnimal), numel(csZone), numel(csPhase));
for iAnimal = 1:numel(csAnimal)
    %----------------------------------------------
    % Plot bars per phase
    subplot(2,5,iAnimal);   
    for iZone = 1:numel(csZone);        
        for iPhase = 1:numel(csPhase)
            eval(sprintf('vsTrialPool = vsTrialPool_%s;', csPhase{iPhase}));
            
            lim = [];
%             if iPhase == 3
%                 lim = [1, 60*2*100]; %limit duration length
%             end
            
            %Load var
            if ~isempty(strfind(strVar, 'RS.'))                
                RS = poolTrials_RS(vsTrialPool, iAnimal-1, lim);                
                vlZ = getZone(RS, iZone);
            elseif ~isempty(strfind(strVar, 'IPI.'))                
                IPI = poolTrials_IPI(vsTrialPool, iAnimal-1, lim);                
                vlZ = getZone(IPI, iZone);
            else
                error('%s not found', strVar);
            end
            eval(sprintf('vrZ = %s;', strVar));
            
            % Zone filtering
            if min(size(vrZ)) > 1
                vrZ1 = vrZ(:,1);
                vrZ2 = vrZ(:,2);
                vrZ = [vrZ1(vlZ), vrZ2(vlZ)];
            elseif numel(vrZ) == numel(vlZ)
                vrZ = vrZ(vlZ);
            else
                disp('not edited');
            end
            
            % compute stats
            mrMean(iZone, iPhase) = fun1(vrZ);
            mrSem(iZone, iPhase) = fun2(vrZ);
            cvZ{iAnimal, iZone, iPhase} = vrZ;
        end        
    end
    plotBarError(mrMean, mrSem, csPhaseColor, csZone);
    AX(end+1) = gca;    
    if iAnimal > 1
        set(gca, 'YTick', []);
    else
        mrMean0 = mrMean;
    end
    title(csAnimal{iAnimal});
    bar(mrMean0, 1, 'FaceColor', 'none', 'EdgeColor', 'k');
    
    %----------------------------------------------
    % Plot differential
    subplot(2,5,iAnimal + 5);
    mrMeanA = [mrMean(:,2) - mrMean(:,1), mrMean(:,3) - mrMean(:,2)];
    mrMeanB = [mrMean(:,2) + mrMean(:,1), mrMean(:,3) + mrMean(:,2)]/2;
    mrMean = mrMeanA./mrMeanB * 100;
    plotBarError(mrMean, [], csDiffColor, csZone);
    AX1(end+1) = gca;    
    if iAnimal > 1
        set(gca, 'YTick', []);
    else
        mrMean1 = mrMean;
    end
    bar(mrMean1, 1, 'FaceColor', 'none', 'EdgeColor', 'k');

end %for

linkaxes(AX);
linkaxes(AX1);
ylabel(AX1(1), '% change');
ylabel(AX(1), strY);

set(AX(1), 'XLim', [1.5 3.5]);
set(AX1(1), 'XLim', [1.5 3.5]);

end %func


function [vlZone, strZone] = getZone(S, iZone)
vlNF = S.vrDf < 14 & S.vrDf >= 3;
vlLM = S.vrD1 < 3 | S.vrD2 < 3 | S.vrD3 < 3 | S.vrD4 < 3;
vlF = S.vrDf < 3;
switch (iZone)
    case 1 %all
%         vlZone = S.tvlZone & ~vlNF & ~vlLM & ~vlF;
%         vlZone = S.vlZone & ~vlLM & ~vlF;
        vlZone = S.vlZone;
        strZone = 'AZ';
    case 2 %LM
        vlZone = vlLM;
        strZone = 'LM<3';
    case 3 %Fc<15
        vlZone = vlNF;
        strZone = 'Fc4~15';
    case 4 %F<3
        vlZone = vlF;
        strZone = 'F<3';
end
vlZone = vlZone(:);
end %func