%S140423_compositeBackground

% this is for the random control
load D140422_Controls;

vsTrials = vsTrialPool_WP;

%%
csAnimals = {'A', 'B', 'C', 'D'};
for iAnimal = 1:4
    S = vsTrials(1);
    mrImg = zeros(size(S.img0));
    nImg = 0;
    for iTrial=1:numel(vsTrials)
        S = vsTrials(iTrial);
        if S.iAnimal == iAnimal
            mrImg = mrImg + double(S.img0);            
            nImg = nImg + 1;
        end
    end

    mrImg = uint8(imadjust(mrImg/nImg/255)*255);
    eval(sprintf('img0_WP_%s = mrImg;', csAnimals{iAnimal}));
end

save img0_WP img0_WP_A img0_WP_B img0_WP_C img0_WP_D;

%%
figure; imshow(imadjust(img0_WP_A)); 
title(sprintf('Mean animal A'));

%% min pooliong
csAnimals = {'A', 'B', 'C', 'D'};
for iAnimal = 1:4
    S = vsTrials(1);
    mrImg = ones(size(S.img0), 'uint8')*255;
    nImg = 0;
    for iTrial=1:numel(vsTrials)
        S = vsTrials(iTrial);
        if S.iAnimal == iAnimal
            mrImg = min(mrImg, imadjust(S.img0));            
        end
    end

    eval(sprintf('img0_WP_%s = mrImg;', csAnimals{iAnimal}));
end

save img0_WP img0_WP_A img0_WP_B img0_WP_C img0_WP_D;

%
figure; imshow(imadjust(img0_WP_A)); 
title(sprintf('Mean animal A'));