function S=appendStruct(S, varargin)
% appends the fields to the existing struct

if numel(varargin) == 1 && isstruct(varargin{1})
    S1 = varargin{1};
    csFieldnames = fieldnames(S1);
    for i=1:numel(csFieldnames)
        fieldname = csFieldnames{i};
        S = setfield(S, fieldname, getfield(S1, fieldname));
    end
else
    for i=1:2:numel(varargin)
        S = setfield(S, varargin{i}, varargin{i+1});
    end
end