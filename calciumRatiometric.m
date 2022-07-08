%% ratiometric ca imaging
% first load file lists into structs
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
%%
% %% background subtraction
% allBGVals = zeros(1,2*length(nFiles));
% for i = 1:2*length(nFiles)
%     if i<=length(nFiles)
%         curIm = imread(nFiles{i});
%     else
%         curIm = imread(dFiles{i-numIms});
%     end
% %     [binValues,binEdges] = histcounts(curIm(:),'BinWidth',1);
% %     binCenters = (binEdges(1:end-1)+binEdges(2:end))/2;
% %     [~,modeIdx] = max(binValues);
% %     fitIdx = modeIdx-50:modeIdx+50;
% %     histPoly = polyfit(binCenters(fitIdx),binValues(fitIdx),2);
% %     curBGVal = roots([histPoly(1)*2, histPoly(2)]);
% %     if curBGVal < 2*mean(curIm(:)) && curBGVal > 0
% %         allBGVals(i) = roots([histPoly(1)*2, histPoly(2)]);
% %     else
% %         allBGVals(i) = mean(curIm(:));
% %     end
%     allBGVals(i) = modalIntensity(curIm);
% end
% nBG = mean(allBGVals(1:numIms));
% dBG = mean(allBGVals(numIms+1:end));
% % nBG = 600;
% % dBG = 320;
 
% %% alt bg from single image
% testIdx = length(nFiles);
% nIm1 = imread(nFiles{testIdx});
% figure
% imagesc(nIm1)
% colormap gray
% title('Crop bg region of numerator image')
% [nIm1,rect] = imcrop;
% close(gcf)
% nBG = modalIntensity(nIm1);
% dIm1 = imread(dFiles{testIdx});
% dIm1 = imcrop(dIm1,rect);
% dBG = modalIntensity(dIm1);

% compute ratio
figure
[saveFile, savePath] = uiputfile('*.tiff', 'Choose location to save ratiometric images');
scaleFactor = 2^10;
colormap gray
for i = 1:numIms
    nIm = double(imread(nFiles{i}));
    nLogic = nIm > nBG + 1500;
    nLogic = imfill(nLogic,'holes');
%     nLogic = nIm < nBG + 10000;
%     nLogic = nIm < nBG + 1500;
    dIm = double(imread(dFiles{i}));
%     dLogic = dIm < dBG + 2000;
    curRatioIm = (nIm-nBG) ./ (dIm-dBG);
%     curRatioIm(nLogic | dLogic) = 0;
    curRatioIm(~nLogic) = 0;
    goodRatioLogic = curRatioIm > .5 & curRatioIm < 1.5;
    goodRatioLogic = imfill(goodRatioLogic,'holes');
    curRatioIm(~goodRatioLogic) = 0;

    curConcIm = 3.892*140*(curRatioIm-.552)./(2.1484-curRatioIm);
%     curRatioIm = 6500./ (dIm-dBG);
    imagesc(curRatioIm)
    caxis([0 5])
    colorbar
    drawnow
    curRatioIm = uint16(curRatioIm*scaleFactor);
    curRatioIm(curRatioIm>=2^16) = 2^16-1;
%     imwrite(curRatioIm,fullfile(savePath,sprintf('%s_%4d.tiff',saveFile(1:end-5),i)))
    imwrite(uint16(curConcIm*10),fullfile(savePath,sprintf('%s_%4d.tiff',saveFile(1:end-5),i)))
end

%% Multiple cell analysis
% identify and track individual cells using blue channel
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
    binIm = nIm > nBG + 2000; % blue cell thresh
%     binIm = nIm > nBG + 10000; % blue cell thresh
    binIm = imfill(binIm, 'holes'); 
    connComps = bwconncomp(binIm); 
    areaProps = regionprops(connComps,'Area','Centroid');
    pxIdxList = connComps.PixelIdxList;
    areas = [areaProps.Area];
    keepAreas = areas > 200 & areas < 2000;
%     keepAreas = areas > 500 & areas < 5000;
    areaProps = areaProps(keepAreas);
    areas = areas(keepAreas);
    pxIdxList = pxIdxList(keepAreas);
    if isempty(areaProps)
        continue
    else
        for j = 1:length(areaProps)
            curLoc = areaProps(j).Centroid;
%             selPx = curRatioIm(pxIdxList{j});
            % select pixels from a radius of 20 pixels
            [xMat,yMat] = meshgrid(1:512,1:512);
            selLogic = sqrt((xMat-curLoc(1)).^2 + (yMat-curLoc(2)).^2) <= 20;
            selPx = curRatioIm(selLogic);
            selPx = selPx(selPx > .5);
            orderedPx = sort(selPx,'ascend');
            curMFI = mean(orderedPx(1:end-5));
%             curIm = curRatioIm(max(curLoc(2)-20,1):min(curLoc(2)+20,512),...
%                 max(curLoc(1)-20,1):min(curLoc(1)+20,512));
%             hold off
%             imagesc(curIm)
%             daspect([1 1 1])
%             hold on
% %             logicIm = zeros(size(curIm));
%             logicIm(curIm>0.7) = 1;
%             b = bwboundaries(logicIm);
%             for k = 1:length(b)
%                 plot(b{k}(:,2),b{k}(:,1))
%             end
%             hold off
            
%             if length(selPx) < 1006
%                 curMFI = mean(orderedPx(1:end-5));
%             else
%                 curMFI = mean(orderedPx(end-1000-5:end-5));
%             end
            placed = false;
            if numCells > 0 && i > 1
                for k = 1:length(cellAnalysis)
                    prevLoc = cellAnalysis(k).loc(i-1,:);
                    if norm(curLoc - prevLoc) < 10 % close -> consider same cell
                        cellAnalysis(k).loc(i,:) = curLoc;
                        cellAnalysis(k).MFI(i) = curMFI;
%                         if isnan(curMFI)
%                             error('uh oh')
%                         end
                        placed = true;
                        break
                    end
                end
            end
            if ~placed % then not identified with cell initiated prev
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

%%
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
        plot(140*3.5*(curMFI-.6)./(2.1-curMFI),'-b')
%         plot(140*3.5*(curMFI-.7)./(2.45-curMFI),'-r')
%         plot(140*4.167*(curMFI-.6)./(2.5-curMFI),'-b')
%         plot(140*6.6*(curMFI-.7)./(4.6-curMFI));
%         plot(140*2.54*(curMFI-.67)./(1.7-curMFI))
%         plot(140*3.7*(curMFI-.59)./(2.2-curMFI))
%         plot(smooth(140*7.14*(curMFI-.7)./(5-curMFI),11))
    end
end

%% quick avg
figure
hold on
allCa = cell(1,length(cellAnalysis(1).MFI));
for i = 1:length(cellAnalysis)
    curMFI = cellAnalysis(i).MFI;
    curMFI(curMFI < .7 | curMFI > 1.5) = nan;
    curCa = 140*3.5*(curMFI-.6)./(2.1-curMFI);
    plot(curCa)
    for j = 1:length(curCa)
        if ~isnan(curCa(j))
            allCa{j} = [allCa{j},curCa(j)];
        end
    end
end
avgCa = zeros(1,length(allCa));
for i = 1:length(allCa)
    avgCa(i) = mean(allCa{i});
end
plot(avgCa,'LineWidth',2)

%% new bg
nBGVals = zeros(numIms,1);
dBGVals = zeros(numIms,1);
for i = 1:numIms
    nIm = double(imread(nFiles{i}));
%     nIm(nIm > 1000) = 0;
%     nBGVals(i) = modalIntensity(nIm);
    nBGVals(i) = mean(nIm(nIm < 1000));
    dIm = double(imread(dFiles{i}));
%     dIm(dIm > 1000) = 0;
%     dBGVals(i) = modalIntensity(dIm);
    dBGVals(i) = mean(dIm(dIm < 1000));
end
figure
plot(nBGVals)
hold on
plot(dBGVals)

%% try with diff bg vals
for i = 1:numIms
    nIm = double(imread(nFiles{i}));
    dIm = double(imread(dFiles{i}));
    curRatioIm = (nIm-nBG) ./ (dIm-dBG);
    for j = 1:length(cellAnalysis)
        if isnan(cellAnalysis(j).loc(i,1))
            cellAnalysis(j).altMFI(i) = nan;
        else
            curLoc = cellAnalysis(j).loc(i,:);
            % select pixels from a radius of 20 pixels
            [xMat,yMat] = meshgrid(1:512,1:512);
            selLogic = sqrt((xMat-curLoc(1)).^2 + (yMat-curLoc(2)).^2) <= 20;
            selPx = curRatioIm(selLogic);
            selPx = selPx(selPx > .5);
            orderedPx = sort(selPx,'ascend');
            cellAnalysis(j).altMFI(i) = mean(orderedPx(1:end-5));
        end
    end
end