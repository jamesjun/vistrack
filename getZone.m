function [vlZone, strZone] = getZone(S, iZone)
vlNF = S.vrDf < 14 & S.vrDf >= 3;
vlLM = S.vrD1 < 3 | S.vrD2 < 3 | S.vrD3 < 3 | S.vrD4 < 3;
vlF = S.vrDf < 3;
if isempty(iZone)
    vlZone = [];
    strZone = '[]';
    return;
end

switch (iZone)
    case 0
        vlZone = [];
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
    case 5 %F<15
        vlZone = vlF | vlNF;
        strZone = 'F<15';        
end
vlZone = vlZone(:);
end %func