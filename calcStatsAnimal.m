function 

mlPath_E1_A = false(size(mlPath_E1));
mlPath_E1_A(:, 1:4) = mlPath_E1(:, 1:4);
meanSpeed_A = bootMean(mrSpeed(mlPath_E1_A));