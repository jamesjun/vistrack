try
    csSettings = get(handles.editSettings, 'String');
    for i=1:numel(csSettings)
        eval(csSettings{i});
    end
catch
    eval('settings');
end