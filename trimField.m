function S = trimField(S, lim, varargin)
% append vector to a struct

%do nothing
if isempty(lim)
    return;
end

vi = lim(1):lim(2);

for ivar = 1:numel(varargin)
    strVar = varargin{ivar};    
    if isfield(S, strVar)        
        vr = getfield(S, strVar);
        if min(size(vr)) == 1
            S = setfield(S, strVar, vr(vi));
        else %matrix
            mr = zeros(numel(vi), size(vr,2));
            for iCol=1:size(mr,2)
                mr(:,iCol) = vr(vi,iCol);
            end
            S = setfield(S, strVar, mr);
        end
    else
        fprintf('var %s doesn''t exist\n', strVar); 
    end
end