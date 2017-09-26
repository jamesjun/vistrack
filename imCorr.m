function mrC = imCorr(mrA, mrB)
% dimension of mrA, mrB must the same

nwin = 21;
nwinh = (nwin-1)/2;
h = size(mrA, 1);
w = size(mrA, 2);

mrA1 = zeros([h+nwinh*2, w+nwinh*2], class(mrA));
mrB1 = zeros([h+nwinh*2, w+nwinh*2], class(mrB));
mrA1(nwinh+1:end-nwinh, nwinh+1:end-nwinh) = mrA;
mrB1(nwinh+1:end-nwinh, nwinh+1:end-nwinh) = mrB;

mlMask = false(h+nwinh*2, nwin); %pad
mrA2 = zeros(nwin^2, h*w);
mrB2 = zeros(nwin^2, h*w);
i=1;
mrC = zeros(h, w);
for ic=(1+nwinh):w-nwinh
    mrA1c = double(mrA1(:, ic-nwinh:ic+nwinh));
    mrB1c = double(mrB1(:, ic-nwinh:ic+nwinh));
    for ir=(1+nwinh):h-nwinh    
        mlMask1 = mlMask;
        mlMask1(ir-nwinh:ir+nwinh, :) = 1;

        mrC(ir-1, ic-1) = cov2(mrA1c(mlMask1), mrB1c(mlMask1));
%         mrA2(:,i) = mrA1c(mlMask1);
%         mrB2(:,i) = mrB1c(mlMask1);        
%         i=i+1;

%         if ic==nwin && ir==nwin
%             disp('debug');
%         end
    end
end
mrC(isnan(mrC)) = 1;
% mrC = reshape(corrMat(mrA2, mrB2, 1), [h, w]);
end %func

function c = cov2(a,b)
    c = mean((a-mean(a)) .* (b-mean(b))) / var(a);
end


function [vrCorr, vrCov] = corrMat(M1, M2, dimm) 
if nargin < 3
    dimm = 1;
end
vrCov = mean(M1 .* M2, dimm) - mean(M1,dimm) .* mean(M2,dimm);
vrCorr = vrCov ./ (nanstd(M1,1,dimm) .* nanstd(M2,1,dimm));
vrCorr(isnan(vrCorr)) = 1;
end