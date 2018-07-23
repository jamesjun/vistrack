try
    csSettings = get(handles.editSettings, 'String');
    for i=1:numel(csSettings)
        eval(csSettings{i});
    end
catch
    eval('settings_vistrack');
end

% Default values
if exist('REPLAY_STEP', 'var') ~= 1, REPLAY_STEP = 4; end