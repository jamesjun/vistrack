function [mrX, mrY] = trial_correctRot(S, angRot)

%rotate background and recalculate the center
%[~, ~, angRot1] = getCoordinatesForAnimal(S.dataID(4));
%img0 = imrotate(S.img0, angRot+angRot1, 'nearest', 'crop');
%xyc = round(size(S.img0)/2);
%rotMat = rotz(angRot1); rotMat = rotMat(1:2, 1:2);
%xy0 = (S.xy0 - S.xyc) * rotMat + S.xy0;

%rotate the path
xyc = S.xy0;
mrX = zeros(size(S.mrX));
mrY = zeros(size(S.mrY));
rotMat = rotz(angRot); rotMat = rotMat(1:2, 1:2);
for iPoint = 1:size(S.mrX,2)
    mrXY = [S.mrX(:,iPoint) - xyc(1), S.mrY(:,iPoint) - xyc(2)] * rotMat;
    mrX(:,iPoint) = mrXY(:,1) + xyc(1);
    mrY(:,iPoint) = mrXY(:,2) + xyc(2);
end