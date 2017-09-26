classdef TrialGroup
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        csAnimals = {'A', 'B', 'C', 'D'};
        angRot = -1.1590; %degs
        nGrid = 20;
        nTime = 25;
    end
    properties
        vTrials; %contains cell of vector of Trial objects
        viAnimals;
        img0;
        xy0; %pooled representation
        vrX, vrY;
        mlMask;
    end
    
    methods
        % Constructor
        % vcName: trial name
        % viAnimals: Animals to pool, optional
        function obj = TrialGroup(vcName, viAnimals)
        if nargin == 0, return; end
        if nargin < 2, viAnimals = []; end
        vcName = upper(vcName);

        switch vcName
            case {'RANDOM', 'LANDMARK', 'NONE'}
            sFnames = struct('RANDOM', 'D141026_RandGroup.mat', ...
                'LANDMARK', 'D140324_LandmarkGroup.mat', ...
                'NONE', 'D141026_NoneGroup.mat');
            eval(sprintf('S = load(sFnames.%s);', vcName));
            obj.viAnimals = [];
            obj.vTrials = Trial();
            for iAnimal = 1:numel(TrialGroup.csAnimals)
                vcAnimal = TrialGroup.csAnimals{iAnimal};
                eval(sprintf('cTrials1 = S.vsTrial_%s;', vcAnimal));
                nTrials1 = numel(cTrials1);
                obj.viAnimals = [obj.viAnimals, iAnimal*ones(1, nTrials1)];
                vTrials1(nTrials1) = Trial();
                for iTrial = 1:nTrials1
                    vTrials1(iTrial) = Trial(cTrials1{iTrial});
                end
                obj.vTrials = [obj.vTrials, vTrials1];
            end
            
            case {'E', 'L', 'P'}
            S = load('D140330_Landmark.mat');
            eval(sprintf('vsTrials = S.vsTrialPool_%s;', vcName));
            vTrials(1, numel(vsTrials)) = Trial();
            for iTrial = 1:numel(vsTrials)
                vTrials(iTrial) = Trial(vsTrials(iTrial));
            end
            obj.viAnimals = [vTrials.iAnimal];
            obj.vTrials = vTrials;    
        end
        
        % Filter by animals
        if ~isempty(viAnimals)
            vl = false(size(obj.viAnimals));
            for iAnimal1 = 1:numel(viAnimals)
                vl = vl | (obj.viAnimals == viAnimals(iAnimal1));
            end
            obj.viAnimals = obj.viAnimals(vl);
            obj.vTrials = obj.vTrials(vl);
        end
        
        % get background image and rotate
        obj.img0 = imrotate(imadjust(obj.vTrials(1).img0), ...
            TrialGroup.angRot, 'nearest', 'crop');
            
        %rotate coordinates
        obj.xy0 = obj.vTrials(1).xy0;
        obj.vrX = [obj.vTrials.vrX];
        obj.vrY = [obj.vTrials.vrY];
        rotMat = rotz(TrialGroup.angRot);    
        rotMat = rotMat(1:2, 1:2);
        mrXY = [obj.vrX(:) - obj.xy0(1), obj.vrY(:) - obj.xy0(2)] * rotMat;
        obj.vrX = mrXY(:,1) + obj.xy0(1);
        obj.vrY = mrXY(:,2) + obj.xy0(2);
        
        % mask
        obj.mlMask = getImageMask(obj.img0, [0 60], 'CENTRE');
        
        obj.vcGroupName = sprintf('%d, ', viAnimals);
        obj.vcGroupName = sprintf('', vcName, );
        end %TrialGroup
          
        
        function obj = filter(varargin)
            P = struct(varargin{:});
            if ~isfield(P, 'iAnimal'), P.iAnimal = []; end
            if ~isfield(P, 'vcLocation'), P.vcLocation = 'ActiveZone'; end
                
            if numel(P.iAnimal) == 1
                obj.vTrials = obj.vTrials(obj.viAnimals == P.iAnimal);
                obj.viAnimals = P.iAnimal * ones(size(obj.vTrials));
            elseif numel(P.iAnimal) > 1
                error('TrialGroup-filter-iAnimal must be scalar');
            end
            
        end %filter
 

        function [RGB, mrPlot] = plotMap(obj, vcMode, lim) 
            if nargin < 3, lim = []; end
            if nargin < 2, vcMode = 'visit'; end

            viX = ceil(obj.vrX / obj.nGrid);
            viY = ceil(obj.vrY / obj.nGrid);
            [h, w] = size(obj.img0);
            h = h / obj.nGrid;
            w = w / obj.nGrid;
            mnVisit = zeros(h, w);
            mnTime = zeros(h, w);
            for iy=1:h
                vlY = (viY == iy);
                for ix=1:w
                    viVisit = find(vlY & (viX == ix));        
                    mnTime(iy,ix) = numel(viVisit);
                    nRepeats = sum(diff(viVisit) < obj.nTime); % remove repeated counts
                    mnVisit(iy,ix) = numel(viVisit) - nRepeats;        
                end
            end

            mrTperV = mnTime ./ mnVisit;
            switch lower(vcMode)
                case 'time'
                    mrPlot = mnTime;
                case 'visit'
                    mrPlot = mnVisit;
                case 'time/visit'
                    mrPlot = mrTperV;
            end

            mnVisit1 = imresize(mrPlot, obj.nGrid, 'nearest');
            mnVisit1(~obj.mlMask) = 0;

            if isempty(lim)
                lim = [min(mnVisit1(:)) max(mnVisit1(:))];
            end

            mrVisit = uint8((mnVisit1 - lim(1)) / diff(lim) * 255);
            RGB = rgbmix(obj.img0, mrVisit, obj.mlMask);
            if nargout == 0, imshow(RGB); end
        end %
    end %methods
end

