%% ratiometric ca imaging analysis
% for Fura Red, the "numerator" refers to the signal following excitation
% by blue light and the "denominator" refers to the signal following
% excitation with cyan light
%
%% first load file lists into structs
[~, nPath] = uigetfile('*.tiff', 'Choose first file of image sequence for numerator');
sequence = dir(fullfile(nPath,'*.tiff'));
nFiles = {sequence.name};
for i = 1:length(nFiles)
    nFiles{i} = fullfile(nPath, nFiles{i});
end
[~, dPath] = uigetfile('*.tiff', 'Choose first file of image sequence for denominator');
sequence = dir(fullfile(dPath,'*.tiff'));
dFiles = {sequence.name};
for i = 1:length(dFiles)
    dFiles{i} = fullfile(dPath, dFiles{i});
end
numIms = min([length(nFiles),length(dFiles)]);
if length(dFiles)~=length(nFiles)
    nFiles = nFiles(1:numIms);
    dFiles = dFiles(1:numIms);
end

%% background subtraction
% can measure background from averaging modal intensities for all images in
% current stack. This will not work well for high densities of cells, best
% to measure modal intensity for a region with unlabeled cells
BGMeasure = questdlg('Quantify background intensity from full stacks?','Measure BG?','Yes','No, bypass','Yes');
if strcmp(BGMeasure, 'Yes')
    allBGVals = zeros(1,2*length(nFiles));
    for i = 1:2*length(nFiles)
        if i<=length(nFiles)
            curIm = imread(nFiles{i});
        else
            curIm = imread(dFiles{i-numIms});
        end
        allBGVals(i) = modalIntensity(curIm);
    end
    nBG = mean(allBGVals(1:numIms));
    dBG = mean(allBGVals(numIms+1:end));
else
    if ~exist('nBG','var')
        BGEntries = inputdlg({'Enter numerator image BG intensity (nBG):', 'Enter denominator image BG intensity (dBG)'},'Enter BG');
        nBG = str2double(BGEntries{1});
        dBG = str2double(BGEntries{2});
    end
end
    
nThresh = nBG + 500; % adjust this value depending on the brightness of cells above bg
dThresh = dBG + 500; % not used here because the denominator images change brightness for FR

%% compute ratio
figure
[saveFile, savePath] = uiputfile('*.tiff', 'Choose location to save ratiometric images');
scaleFactor = 2^10; % used to convert ratio into a 16 bit image (e.g. intensity of 1024 corresponds to ratio of 1)
colormap gray
for i = 1:numIms
    nIm = double(imread(nFiles{i}));
    nLogic = nIm > nThresh; 
    nLogic = imfill(nLogic,'holes');
    dIm = double(imread(dFiles{i}));
%     dLogic = dIm > dThresh;
    curRatioIm = (nIm-nBG) ./ (dIm-dBG);
%     curRatioIm(~nLogic | ~dLogic) = 0; % see comment about dThresh above
    curRatioIm(~nLogic) = 0;
    goodRatioLogic = curRatioIm > .5 & curRatioIm < 2;
    goodRatioLogic = imfill(goodRatioLogic,'holes');
    curRatioIm(~goodRatioLogic) = 0; % only consider ratios between .5 and 2

%     curConcIm = 3.892*140*(curRatioIm-.552)./(2.1484-curRatioIm); % can
%     convert to concentrations if desired
    imagesc(curRatioIm)
    caxis([0 5])
    colorbar
    drawnow
    curRatioIm = uint16(curRatioIm*scaleFactor);
    curRatioIm(curRatioIm>=2^16) = 2^16-1;
    imwrite(curRatioIm,fullfile(savePath,sprintf('%s_%4d.tiff',saveFile(1:end-5),i)))
%     imwrite(uint16(curConcIm*10),fullfile(savePath,sprintf('%s_%4d.tiff',saveFile(1:end-5),i)))
end

%% Multiple cell analysis
% identify and track individual cells using blue (numerator) channel
figure
colormap gray
cellAnalysis = struct('loc',[],'MFI',[],'Ca',[]); % for storing data by cell
numCells = 0;
for i = 1:numIms
    nIm = double(imread(nFiles{i}));
    dIm = double(imread(dFiles{i}));
    curRatioIm = (nIm-nBG) ./ (dIm-dBG);
    imagesc(nIm)
    hold on
    binIm = nIm > nBG + 2000; % threshold to identify single cells - should be adjusted according to the brightness of the cells above the background
    binIm = imfill(binIm, 'holes'); 
    connComps = bwconncomp(binIm); 
    areaProps = regionprops(connComps,'Area','Centroid');
    pxIdxList = connComps.PixelIdxList;
    areas = [areaProps.Area];
    keepAreas = areas > 200 & areas < 2000; % this size filter should be adjusted based on the scale value and the size of the cells being considered
%     keepAreas = areas > 500 & areas < 5000;
    areaProps = areaProps(keepAreas);
    areas = areas(keepAreas);
    pxIdxList = pxIdxList(keepAreas);
    if isempty(areaProps) % no cells identified
        continue
    else
        for j = 1:length(areaProps)
            curLoc = areaProps(j).Centroid;
            % select pixels from a radius of 20 pixels
            [xMat,yMat] = meshgrid(1:512,1:512);
            selLogic = sqrt((xMat-curLoc(1)).^2 + (yMat-curLoc(2)).^2) <= 20;
            selPx = curRatioIm(selLogic);
            selPx = selPx(selPx > .5);
            orderedPx = sort(selPx,'ascend');
            curMFI = mean(orderedPx(1:end-5));
            % could optionally only use a set number of pixels within this
            % region (1000 in the commented-out code below)
%             if length(selPx) < 1006
%                 curMFI = mean(orderedPx(1:end-5));
%             else
%                 curMFI = mean(orderedPx(end-1000-5:end-5));
%             end

            % now, determine whether the current cell corresponds to one
            % previously identified (here, if a centroid is closer than 10
            % pixels to a cell body in the previous frame, they are
            % identified as the same cell)
            placed = false;
            if numCells > 0 && i > 1
                for k = 1:length(cellAnalysis)
                    prevLoc = cellAnalysis(k).loc(i-1,:);
                    if norm(curLoc - prevLoc) < 10 % close -> consider same cell
                        cellAnalysis(k).loc(i,:) = curLoc;
                        cellAnalysis(k).MFI(i) = curMFI;
                        placed = true;
                        break
                    end
                end
            end
            if ~placed % then not identified with cell initiated prev, start a new cell in the cellAnalysis structure
                numCells = numCells + 1;
                cellAnalysis(numCells).loc = nan(numIms,2);
                cellAnalysis(numCells).MFI = nan(numIms,1);
                cellAnalysis(numCells).loc(i,:) = curLoc;
                cellAnalysis(numCells).MFI(i,:) = curMFI;
            end
        end
    end
    for k = 1:length(areaProps)
        curCentroid = areaProps(k).Centroid;
        plot(curCentroid(1),curCentroid(2),'r*')
    end
    title(sprintf('Cells at frame %d',i))
    drawnow
    hold off
end

%% optionally, delete any cell analyses shorter than 15 frames
keepLogic = true(numCells,1);
for i = 1:numCells
    if sum(~isnan(cellAnalysis(i).MFI)) < 15
        keepLogic(i) = false;
    end
end
cellAnalysis = cellAnalysis(keepLogic);

%% make plots
figure
hold on
for i = 1:length(cellAnalysis)
    curMFI = cellAnalysis(i).MFI;
%     curMFI = curMFI(~isnan(curMFI));
    if any(curMFI > .7)
%         plot(curMFI)
        plot(3.892*140*(curMFI-.552)./(2.1484-curMFI),'-b') % includes calibrated values from PLB-985 calibration (see EAF dissertation)
    end
end