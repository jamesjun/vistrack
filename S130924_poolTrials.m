%

dname = 'C:\Expr\2013.2S_CycleW_probe\';
csFname_A = {'W06A4p_Track', 'W07A3p_Track', 'W08A2p_Track', 'W09A1p_Track'};
csFname_B = {'W06B4p_Track', 'W07B3p_Track', 'W08B2p_Track', 'W09B1p_Track'};
csFname_C = {'W06C4p_Track', 'W07C3p_Track', 'W08C2p_Track'};
csFname_D = {'W06D4p_Track', 'W07D3p_Track', 'W08D2p_Track', 'W09D1p_Track'};

%%
fishID = {'A', 'B', 'C', 'D'};
for iFish=1:4
    eval(sprintf('csFname = csFname_%s;', fishID{iFish}));
    cmVISITCNT = {};  cmTIMECNT = {};
    for i=1:numel(csFname)
        handles = load([dname csFname{i}]);
        [cmVISITCNT{i}, cmTIMECNT{i}] = calcVisitDensity_aux(handles);
    end
    eval(sprintf('cmVISITCNT_%s = cmVISITCNT;', fishID{iFish}));
    eval(sprintf('cmTIMECNT_%s = cmTIMECNT;', fishID{iFish}));
end

%
save D130924_CycleW cmTIMECNT_* cmVISITCNT_*;

%% load images
handles = load([dname, csFname_A{1}]);
figure; imshow(handles.img0);

%%

%%
handles = load([dname, csFname_B{3}]);
I0 = imresize(handles.img0, .25);

figure; 
subplot 122; 
imshow(rgbmix(I0, imgray2rgb((cmVISITCNT_B{3}))));   title('Time density map');
% 
subplot 121; 
imshow(handles.img0);