function [ S ] = makeStruct( varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name

S=[];

for i=1:nargin
    S = setfield(S, inputname(i), varargin{i});
end

end