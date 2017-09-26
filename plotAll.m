function plotAll(csTrials, csCmd, viZone, csX, nCols)
% csCmd: pair: command, ylabel
% iZone: optional. default 1
% csX: optional. default: {E,L,P}
%      pair: condition, XTickLabel

vcAnimal = 'o^sd'; %animal's shape
vcPhase = 'rbg';
csLine = ':';
mrMean = zeros(4,3);
mrSem = zeros(4,3);

if nargin < 3
    viZone = 1;
    strZone = '';
end
if isempty(viZone)
    viZone = 1;
end
if nargin < 4
    csX = [];
end
if nargin < 5
    nCols = [];
end
nPlot = size(csCmd,1);
if numel(viZone) ~= nPlot
    viZone = ones(nPlot,1) * viZone(1);
end

vlDep = @(S)(S.vrD1 <= 3 & differentiate3(S.vrD1) > 0) |...
               (S.vrD2 <= 3 & differentiate3(S.vrD2) > 0) |...
               (S.vrD3 <= 3 & differentiate3(S.vrD3) > 0) |...
               (S.vrD4 <= 3 & differentiate3(S.vrD4) > 0);
         
vlApp = @(S)(S.vrD1 <= 3 & differentiate3(S.vrD1) < 0) |...
               (S.vrD2 <= 3 & differentiate3(S.vrD2) < 0) |...
               (S.vrD3 <= 3 & differentiate3(S.vrD3) < 0) |...
               (S.vrD4 <= 3 & differentiate3(S.vrD4) < 0);
nansem = @(x)nanstd(x) / sqrt(sum(~isnan(x)));
           
nPlot = size(csCmd,1);
switch size(csCmd,2)
    case 2
        fun1 = @(x)nanmean(x);
        fun2 = @(x)nansem(x);
    case 3
        fun1 = [];
        fun2 = @(x)nansem(x);  
end
if isempty(nCols)
    nCols = size(csCmd, 1);
end
nRows = ceil(nPlot/nCols);
nAnimals = numel(vcAnimal);
if isempty(csX)
    nX = numel(csTrials);
    csXstr = {'E', 'L', 'P'};
elseif min(size(csX)) == 1
    csXstr = csX;
    nX = numel(csX);
    csX = [];
else
    nX = size(csX, 1);    
    csXstr = csX(:,2);    
end
vrX = 1:nX;
figure;

for iCmd=1:nPlot
    subplot(nRows, nCols, iCmd); hold on;
    strCmd = csCmd{iCmd,1};
    if size(csCmd,2) >= 3
        eval(sprintf('fun1 = %s;', csCmd{iCmd,3}));
    end
    if size(csCmd,2) >= 4
        eval(sprintf('fun2 = %s;', csCmd{iCmd,4}));
    end
    mrY = zeros(nAnimals, nX);
    mrE = zeros(nAnimals, nX);
    for iPhase = 1:numel(csTrials)
        vsTrialPool = csTrials{iPhase};
        for iAnimal = 1:nAnimals               
            if ~isempty(strfind(strCmd, 'RS'))                
                RS = poolTrials_RS(vsTrialPool, iAnimal, [], viZone(iCmd));                
                [vlZ0, strZone] = getZone(RS, viZone(iCmd));
            elseif ~isempty(strfind(strCmd, 'IPI'))                
                IPI = poolTrials_IPI(vsTrialPool, iAnimal, [], viZone(iCmd));                
                [vlZ0, strZone] = getZone(IPI, viZone(iCmd));
            else
                error('%s not found', strCmd);
            end
            eval(sprintf('vrZ = %s;', strCmd));            
            if ~isempty(csX)
                for iX=1:nX
                    eval(sprintf('vlZ = vlZ0 & %s;', csX{iX,1}));
                    [vrY(iX), vrE(iX)] = calcZstats(vrZ, vlZ, fun1, fun2);
                end
            else
                [mrY(iAnimal,iPhase), mrE(iAnimal,iPhase)] = ...
                    calcZstats(vrZ, vlZ0, fun1, fun2);
            end             
            %plot here
            if ~isempty(csX)
                errorbar(vrX, vrY, vrE, [vcPhase(iPhase), vcAnimal(iAnimal), csLine]);
            end            
        end %iAnimal                 
    end %iPhase
    if isempty(csX)
        for iAnimal=1:nAnimals
            errorbar(vrX, mrY(iAnimal,:), mrE(iAnimal,:), ...
                    ['k', vcAnimal(iAnimal), csLine]);
        end
    end
              
    %plot
    ylabel(csCmd{iCmd,2});
    set(gca, 'XTick', 1:nX);
    set(gca, 'XLim', [.5, nX+.5]);
    set(gca, 'XTickLabel', csXstr);        
end %for

suptitle(strZone);

end %func


function [v1, v2] = calcZstats(vrZ, vlZ, fun1, fun2)
% Zone filtering
if ~isempty(vlZ)
    if min(size(vrZ)) > 1
        vrZ1 = vrZ(:,1);
        vrZ2 = vrZ(:,2);
        vrZ = [vrZ1(vlZ), vrZ2(vlZ)];
    elseif numel(vrZ) == numel(vlZ)
        vrZ = vrZ(vlZ);
    else
        warning('vlZ not edited');
    end
end

% compute stats
if nargin < 3 %MEAN
    v1 = nanmean(vrZ);
else
    v1 = fun1(vrZ);
end
if nargin < 4 %SEM
    v2 = nanstd(vrZ) / sqrt(sum(~isnan(vrZ)));
else
    v2 = fun2(vrZ);
end
end