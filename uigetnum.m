function vrNum = uigetnum(csPrompt, def_val)
% Examples
% vrNum = uigetnum('num1');
% vrNum = uigetnum('num1', 1);
% vrNum = uigetnum({'num1', 'num2'}, [1,2])

if nargin<1, csPrompt = {}; end
if nargin<2, def_val = ''; end
if ischar(csPrompt), csPrompt = {csPrompt}; end
if ischar(def_val)
    cs_def_val = {def_val};
elseif iscell(def_val)
    cs_def_val = def_val;
else
    cs_def_val = arrayfun(@(x)num2str(x), def_val, 'UniformOutput', 0);
end

csAns = inputdlg(csPrompt, '', 1, cs_def_val);
if isempty(csAns), vrNum=[]; return; end

vrNum = nan(size(csPrompt));
for iAns = 1:numel(csAns)
    try
        vrNum(iAns) = str2num(csAns{iAns});
    catch
        ;
    end
end
end %func