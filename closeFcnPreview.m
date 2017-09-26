function closeFcnPreview(hObject, event)

S = get(hObject, 'UserData');
try
    stop(S);
    delete(S);
catch
    disp(lasterr);
end
delete(hObject);