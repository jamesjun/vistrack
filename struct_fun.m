% 7/20/2018
% James Jun
function varargout = struct_fun(varargin)
% S_save = struct_copy_(handles, csField)


if nargin==0
    vcCmd = 'help'; 
else
    vcCmd = varargin{1};
end

switch vcCmd
    case 'help', help_(); 
    case 'copy', copy_();
    case 'get', get_();
    case 'set', set_();
end %switch

end %func


