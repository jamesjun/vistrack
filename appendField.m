function S = appendField(S, varargin)
% append vector to a struct

try
for ivar = 1:numel(varargin)
    strVar = inputname(ivar+1);
    vr = varargin{ivar};
    vr = vr(:);
    if isfield(S, strVar)        
        vr = [getfield(S, strVar); vr];        
    end
    S = setfield(S, strVar, vr);
end

catch
   disp('blah'); 
end