%% compile ratio data, align
% figure
[MFIFile,MFIPath] = uigetfile('*.mat','Choose .mat file with MFI data');
MFIData = load(fullfile(MFIPath,MFIFile));
MFIData = MFIData.cellAnalysis;
[timeFile, timePath] = uigetfile('*.txt','Choose fluorescence time .txt file');
timeDataRaw = load(fullfile(timePath, timeFile));
numFrames = length(MFIData(1).MFI);
timeData = (timeDataRaw(1:numFrames,1)-timeDataRaw(1,1)) / 1000; % convert to rel time, seconds
timeInterp = -50:2:200;
sortedMFI = cell(1,length(timeInterp));
MFIStored = cell(1,length(MFIData));
allMFIFig = figure;
curSite = '';
date = inputdlg('Enter date:');
date = date{1};
[~, imPath] = uigetfile('*.tiff', 'Choose first file of image sequence for sample seq');
sequence = dir(fullfile(imPath,'*.tiff'));
imFiles = {sequence.name};
for i = 1:length(imFiles)
    imFiles{i} = fullfile(imPath, imFiles{i});
end
imSize = size(imread(imFiles{1}));

ratiometric = false;
if ~ratiometric
    [BGFile,BGPath] = uigetfile('*.mat','Choose .mat file with BG value');
    BGLoad = load(fullfile(BGPath,BGFile));
    areaAtFlux = nan(length(MFIData),1);
    [timeFile, timePath] = uigetfile('*.txt','Choose RICM time .txt file');
    timeDataRICM = load(fullfile(timePath, timeFile));
    timeData = (timeDataRaw(1:numFrames,1)-timeDataRICM(1,1)) / 1000; % convert to rel time, seconds
end

for i = 1:length(MFIData)
    MFIFig = figure;
    curMFI = MFIData(i).MFI;
    selPts = ~isnan(curMFI);
    curMFI = curMFI(selPts);
    BG = 0;
    if ratiometric
%         curMFI = 3.5*140*(curMFI-.6)./(2.1-curMFI); % convert to Ca
        curMFI = 3.892*140*(curMFI-.552)./(2.1484-curMFI);
%     curMFI = 3.5*140*(curMFI-.7)./(2.45-curMFI); % convert to Ca
    else
        curMFI = curMFI - BGLoad.BG;
    end
    curTime = timeData(selPts);
    [curMFI_noOL,logicOL] = rmoutliers(curMFI,'movmed',9);
    curTime = curTime(~logicOL);
    curMFI = curMFI_noOL;
%     curMFI = interp1(curTime(~logicOL),curMFI_noOL,curTime);
    
    
    plot(curTime,curMFI)
    xlabel('rel time (s)')
    ylabel('MFI (au)')
    title(sprintf('Cell #%d',i))
    hold on
    isFlux = questdlg('Is there a flux?','Flux?','Yes','No','Yes');
    if strcmp(isFlux,'No')
        close(gcf)
        continue
    end
  
    locVals = MFIData(i).loc;
    figure
    colormap gray
    selIdx = find(selPts);
    imStack = zeros(51,51,length(selIdx));
    for j = selIdx(1):selIdx(end)
        curIm = imread(imFiles{j});
        curLoc = round(locVals(j,:));
        rowIdx = curLoc(2)-25:curLoc(2)+25;
        rowIdx = rowIdx(rowIdx > 0 & rowIdx <= imSize(1));
        colIdx = curLoc(1)-25:curLoc(1)+25;
        colIdx = colIdx(colIdx > 0 & colIdx <= imSize(2));
        croppedIm = zeros(51,51);
        localRow = rowIdx - curLoc(2)+26;
        localCol = colIdx - curLoc(1)+26;
        croppedIm(localRow,localCol) = curIm(rowIdx,colIdx);
        imStack(:,:,selIdx==j) = croppedIm/50000;
        
        if j == selIdx(1) && ~ratiometric
            figure
%             imagesc(curIm)
            imagesc(imrotate(curIm,180))
            colormap gray
            hold on
%             patch('XData',[colIdx(1), colIdx(end), colIdx(end), colIdx(1)],...
%                 'YData',[rowIdx(1), rowIdx(1), rowIdx(end), rowIdx(end)],...
%                 'EdgeColor','y','FaceColor','none','LineWidth',2);
            patch('XData',[size(curIm,2)-colIdx(1)+1, size(curIm,2)-colIdx(end)+1, size(curIm,2)-colIdx(end)+1, size(curIm,2)-colIdx(1)+1],...
                'YData',[size(curIm,1)-rowIdx(1)+1, size(curIm,1)-rowIdx(1)+1, size(curIm,1)-rowIdx(end)+1, size(curIm,1)-rowIdx(end)+1],...
                'EdgeColor','y','FaceColor','none','LineWidth',2);
            areaDataLogic = questdlg('Is there area data for this cell?','Area data?','Yes','No','Yes');
            if strcmp(areaDataLogic, 'Yes')
                [areaFile,areaPath] = uigetfile('*.txt','Load area file');
                areaData = load(fullfile(areaPath,areaFile));
            else
                areaData = [curTime, nan(size(curTime))];
            end
            close(gcf)
        end
        
        imagesc(croppedIm)
        drawnow
    end
    doesSpread = questdlg('Does this cell spread?','Spread?','Yes','No','Yes');
    if strcmp(doesSpread,'No') && strcmp(isFlux,'No')
        close(gcf)
        close(MFIFig)
        continue
    elseif strcmp(doesSpread,'No')
        close(gcf)
        startSpreadTime = NaN;
    elseif strcmp(doesSpread,'Yes') && strcmp(isFlux,'No')
        save('imStack.mat','imStack')
        waitfor(chooseSpreadingFrame)
        startSpreading = load('spreadingFrame.mat');
        startSpreading = startSpreading.spreadingFrame;
        delete('imStack.mat')
        delete('spreadingFrame.mat')
        close(gcf)
        startSpreadTime = curTime(startSpreading);
    else
        close(gcf)
        startSpreadTime = 1;
    end
    
%     isFlux = questdlg('Is there a flux?','Flux?','Yes','No','Yes');
%     if strcmp(isFlux,'No')
%         close(gcf)
%         continue
%     end

    curCode = sprintf('%d',i);
    if strcmp(isFlux,'Yes')
        uiwait(msgbox('Select start and end of baseline, then start and end of peak'))
        [BaseTime1,~] = selectPoint(curTime,curMFI);
        [BaseTime2,~] = selectPoint(curTime,curMFI);
        startBase = find(abs(BaseTime1-curTime)<1e-12);
        endBase = find(abs(BaseTime2-curTime)<1e-12);
        BaseMFI = mean(curMFI(startBase:endBase));
        plot([BaseTime1,BaseTime2],[BaseMFI,BaseMFI],'-k')
        [PeakTime1,~] = selectPoint(curTime,curMFI);
        [PeakTime2,~] = selectPoint(curTime,curMFI);
        startPeak = find(abs(PeakTime1-curTime)<1e-12);
        endPeak = find(abs(PeakTime2-curTime)<1e-12);
        if length(startPeak:endPeak) >= 3
            peakFit = polyfit(curTime(startPeak:endPeak),curMFI(startPeak:endPeak),2);
            dMFI = polyder(peakFit);
            peakTime = roots(dMFI);
            if peakTime > PeakTime1 && peakTime < PeakTime2 && peakFit(1) < 0
                PeakMFI = polyval(peakFit,peakTime);
                plot(curTime(startPeak:endPeak),polyval(peakFit,curTime(startPeak:endPeak)),'-k')
            else
                PeakMFI = mean(curMFI(startPeak:endPeak));
                plot([PeakTime1,PeakTime2],[PeakMFI,PeakMFI],'-k')
            end
        else
            PeakMFI = mean(curMFI(startPeak:endPeak));
            plot([PeakTime1,PeakTime2],[PeakMFI,PeakMFI],'-k')
        end
        baseline = BaseMFI;
        fluxRatio = (PeakMFI-BG) / (BaseMFI-BG);
        baseCa = BaseMFI - BG;
        peakCa = PeakMFI - BG;
        % find calcium burst timing
        [startTime, peakTimes, automatic, extraVals] = findFluxPoint3(curTime,curMFI,BG,baseline);
        %     BG = extraVals(1);
        %     fluxRatio = extraVals(2);
        if length(peakTimes) > 5
            warning('more than 5 peaks')
            peakTimes = peakTimes(1:5);
        else
            catVals = zeros(1,5-length(peakTimes));
            peakTimes = [peakTimes,catVals];
        end
        outputCell(1:2) = {date, curCode};
        outputCell(3:10) = {startTime, peakTimes(1),peakTimes(2),peakTimes(3),peakTimes(4),peakTimes(5),fluxRatio,BG};
        outputCell(11:13) = {startSpreadTime, baseCa, peakCa};
    else
        outputCell(1:2) = {date, curCode};
        outputCell(3:10) = {NaN, NaN, NaN, NaN, NaN, NaN, 1, BG};
        outputCell(11:13) = {startSpreadTime, baseCa, peakCa};
        startTime = startSpreadTime;
    end
    if ~ratiometric
        if isnan(areaData(1,2))
            outputCell{1} = sprintf('%s - MFI only', date);
        else
            outputCell{1} = sprintf('%s - %s', date, areaFile(1:end-4));
        end
        if startTime < areaData(1,1)
            startRad = sqrt(areaData(1,2)/pi);
            if (areaData(1,1)-startTime) < startRad/.08
                outputCell{10} = pi*(startRad - (areaData(1,1)-startTime)*.08)^2;
            else
                outputCell{10} = 0;
            end
        else
            outputCell{10} = interp1(areaData(:,1),areaData(:,2),startTime);
        end
    end
    choice = questdlg('Save new file with header or append to old?', 'New file?', 'Append', 'New file', 'Append');
    switch choice
        case 'Append'
            [file,path] = uigetfile('*.txt', 'Choose existing data matrix file');
            fid = fopen(fullfile(path,file),'a');
            fprintf(fid, '%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', outputCell{1:13});
            fclose(fid);
        case 'New file'
            [file,path] = uiputfile('*.txt', 'Choose location to save new data matrix file');
            fid = fopen(fullfile(path,file),'w');
            header(1:2) = {'Date', 'Cell'};
            header(3:13) = {'Start Flux Time', 'Peak Time 1', 'Peak Time 2','Peak Time 3', 'Peak Time 4', 'Peak Time 5',...
                'Flux Ratio','BG MFI', 'Start Spreading Time','Resting Calcium','Peak Calcium'};
            if ~ratiometric
                header{10} = 'Area at Flux';
            end
            fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n', header{1:13});
            fprintf(fid, '%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', outputCell{1:13});
            fclose(fid);
    end
    close(gcf)
    if ~ratiometric
        curMFI = curMFI / baseCa;
    end
    curTime = curTime-startTime;
%     indexGap = round(50 / mean(diff(curTime)));
%     startBurstIdx = find(abs(curTime)<1e-3);
%     baselineIdx = max(1,startBurstIdx-indexGap-1):startBurstIdx;
    figure(allMFIFig)
    plot(curTime,curMFI)
    hold on
    MFIStored{i} = [curTime, curMFI];
    startInterpTime = max([ceil(curTime(1)/2)*2,-50]);
    endInterpTime = min([floor(curTime(end)/2)*2,200]);
    curTimeInterp = startInterpTime:2:endInterpTime;
    MFIInterp = interp1q(curTime,curMFI,curTimeInterp');
    for j = 1:length(MFIInterp)
        [~,curIdx] = min(abs(timeInterp-curTimeInterp(j)));
        sortedMFI{curIdx} = [sortedMFI{curIdx},MFIInterp(j)];
    end
end
xlim([-50 200])
meanMFI = zeros(1,length(sortedMFI));
stdMFI = zeros(1,length(sortedMFI));
semMFI = zeros(1,length(sortedMFI));
for i = 1:length(sortedMFI)
    meanMFI(i) = mean(sortedMFI{i});
    stdMFI(i) = std(sortedMFI{i});
    semMFI(i) = std(sortedMFI{i}) / sqrt(length(sortedMFI{i}));
end
plot(timeInterp,meanMFI,'LineWidth',2)



%% load saved data analysis
[MFIFile,MFIPath] = uigetfile('*.mat','Choose .mat file with MFI data');
MFIData = load(fullfile(MFIPath,MFIFile));
MFIData = MFIData.cellAnalysis;
[timeFile, timePath] = uigetfile('*.txt','Choose time .txt file');
timeData = load(fullfile(timePath, timeFile));
numFrames = length(MFIData(1).MFI);
timeData = (timeData(1:numFrames,1)-timeData(1,1)) / 1000; % convert to rel time, seconds
timeInterp = -50:2:200;
[alignFile, alignPath] = uigetfile('*.txt','Choose alignment analysis .txt file');
alignData = readcell(fullfile(alignPath,alignFile));
numCells = size(alignData,1) - 1; % first row is header
sortedMFI = cell(1,length(timeInterp));
MFIStored = cell(1,length(MFIData));
figure
hold on

[~, imPath] = uigetfile('*.tiff', 'Choose first file of image sequence for sample seq');
sequence = dir(fullfile(imPath,'*.tiff'));
imFiles = {sequence.name};
for i = 1:length(imFiles)
    imFiles{i} = fullfile(imPath, imFiles{i});
end
imSize = size(imread(imFiles{1}));
startSpreadingVec = zeros(numCells,1);

for i = 1:numCells
    curCellNum = alignData{i+1,2};
    curMFI = MFIData(curCellNum).MFI;
    selPts = ~isnan(curMFI);
    curMFI = curMFI(selPts);
%     curMFI = 3.5*140*(curMFI-.6)./(2.1-curMFI); % convert to Ca
    curMFI = 3.892*140*(curMFI-.552)./(2.1484-curMFI);
    curTime = timeData(selPts);
    selIdx = find(selPts);
    
    % opt find start spreading
%     locVals = MFIData(curCellNum).loc;
%     imStack = zeros(51,51,length(selIdx));
%     figure
%     set(gcf,'Position',[400 200 1000 500])
%     subplot(1,2,1)
%     plot(curTime,curMFI)
%     title(sprintf('Cell %d MFI',curCellNum))
%     subplot(1,2,2)
%     colormap gray
%     for j = 1:length(selIdx)
%         jIdx = selIdx(j);
%         curIm = imread(imFiles{jIdx});
%         curLoc = round(locVals(jIdx,:));
%         rowIdx = curLoc(2)-25:curLoc(2)+25;
%         rowIdx = rowIdx(rowIdx > 0 & rowIdx <= imSize(1));
%         colIdx = curLoc(1)-25:curLoc(1)+25;
%         colIdx = colIdx(colIdx > 0 & colIdx <= imSize(2));
%         croppedIm = zeros(51,51);
%         localRow = rowIdx - curLoc(2)+26;
%         localCol = colIdx - curLoc(1)+26;
%         croppedIm(localRow,localCol) = curIm(rowIdx,colIdx);
%         imStack(:,:,selIdx==jIdx) = croppedIm/50000;
%         imagesc(croppedIm)
%         drawnow
%     end
%     doesSpread = questdlg('Does this cell spread?','Spread?','Yes','No','Yes');
%     if strcmp(doesSpread,'No')
%         alignData{i+1,11} = NaN;
%     else
%         save('imStack.mat','imStack')
%         waitfor(chooseSpreadingFrame)
%         startSpreading = load('spreadingFrame.mat');
%         startSpreading = startSpreading.spreadingFrame;
%         delete('imStack.mat')
%         delete('spreadingFrame.mat')
%         startSpreadingVec(i) = curTime(startSpreading);
%         alignData{i+1,11} = startSpreadingVec(i);
%     end
%     close(gcf)
    
    if isnan(alignData{i+1,3})
        curTime = curTime-alignData{i+1,11};
    else
        curTime = curTime-alignData{i+1,3};
    end
%     curTime = curTime-alignData{i+1,11};
    plot(curTime,curMFI)
    MFIStored{i} = [curTime, curMFI];
    startInterpTime = max([ceil(curTime(1)/2)*2,-50]);
    endInterpTime = min([floor(curTime(end)/2)*2,200]);
    curTimeInterp = startInterpTime:2:endInterpTime;
    MFIInterp = interp1q(curTime,curMFI,curTimeInterp');
    for j = 1:length(MFIInterp)
        [~,curIdx] = min(abs(timeInterp-curTimeInterp(j)));
        sortedMFI{curIdx} = [sortedMFI{curIdx},MFIInterp(j)];
    end
end
meanMFI = zeros(1,length(sortedMFI));
stdMFI = zeros(1,length(sortedMFI));
semMFI = zeros(1,length(sortedMFI));
for i = 1:length(sortedMFI)
    meanMFI(i) = mean(sortedMFI{i});
    stdMFI(i) = std(sortedMFI{i});
    semMFI(i) = std(sortedMFI{i}) / sqrt(length(sortedMFI{i}));
end
plot(timeInterp,meanMFI,'LineWidth',2)

%% save new file after adding spreading time
[file,path] = uiputfile('*.txt', 'Choose location to save new data matrix file');
fid = fopen(fullfile(path,file),'w');
header(1:2) = {'Date', 'Cell'};
header(3:13) = {'Start Flux Time', 'Peak Time 1', 'Peak Time 2','Peak Time 3', 'Peak Time 4', 'Peak Time 5',...
    'Flux Ratio','BG MFI', 'Start Spreading Time','Resting Calcium','Peak Calcium'};
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n', header{1:13});
for i = 1:numCells
    fprintf(fid, '%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', alignData{i+1,1:13});
end
fclose(fid);

%% compile condition
MFIDataSummary = struct;
timeInterp = -50:2:200;
sortedMFI = cell(1,length(timeInterp));
figure
hold on
baselines = [];
peaks = [];

loadMore = true;
while loadMore
    [MFIFile,MFIPath] = uigetfile('*.mat','Choose .mat file with MFI data');
    MFIData = load(fullfile(MFIPath,MFIFile));
    MFIData = MFIData.cellAnalysis;
    [timeFile, timePath] = uigetfile(strcat(MFIPath,'*.txt'),'Choose time .txt file');
    timeData = load(fullfile(timePath, timeFile));
    numFrames = length(MFIData(1).MFI);
    timeData = (timeData(1:numFrames,1)-timeData(1,1)) / 1000; % convert to rel time, seconds
    [alignFile, alignPath] = uigetfile(strcat(MFIPath,'*.txt'),'Choose alignment analysis .txt file');
    alignData = readcell(fullfile(alignPath,alignFile));
    numCells = size(alignData,1) - 1; % first row is header
    MFIStored = cell(1,length(MFIData));
    for j = 1:numCells
        curDate = alignData{j+1,1};
        curCellNum = alignData{j+1,2};
        curMFI = MFIData(curCellNum).MFI;
%         curMFI = MFIData(curCellNum).altMFI';
        selPts = ~isnan(curMFI);
        curMFI = curMFI(selPts);
%         curMFI = 3.8*140*(curMFI-.55)./(2.1-curMFI); % convert to Ca
%         curMFI = 3.5*140*(curMFI-.6)./(2.1-curMFI); % convert to Ca
          curMFI = 3.892*140*(curMFI-.552)./(2.1484-curMFI);
%         curMFI = smooth(curMFI,'rlowess');
        curTime = timeData(selPts);
        if isnan(alignData{j+1,3})
            zeroTime = alignData{j+1,11};
            curTime = curTime-zeroTime;
        else
            zeroTime = alignData{j+1,3};
            curTime = curTime-zeroTime;
        end
        %     curTime = curTime-alignData{j+1,11};
        plot(curTime,curMFI)
        % get resting and peak
        
        startFlux = 0; % start flux or start spreading, depending on cell
        [~,fluxIdx] = min(abs(curTime - startFlux));
        
        if size(alignData,2) <= 11
            %         startBase = max([1,fluxIdx-5]);
            %         smoothMFI = smooth(curMFI,'rlowess');
            %         baselines = vertcat(baselines, mean(smoothMFI(startBase:fluxIdx)));
            %         diffSmoothMFI = diff(smoothMFI);
            %         peakIdx = find(diffSmoothMFI(fluxIdx:end)<0,1,'first');
            %         peakIdx = peakIdx + fluxIdx-1;
            %         if isempty(peakIdx)
            %             warning('no peak')
            %             peaks = vertcat(peaks, nan);
            %         end
            %         peaks = vertcat(peaks, curMFI(peakIdx));
            figure
            plot(curTime,curMFI)
            hold on
            uiwait(msgbox('Select start and end of baseline, then start and end of peak'))
            [BaseTime1,~] = selectPoint(curTime,curMFI);
            [BaseTime2,~] = selectPoint(curTime,curMFI);
            startBase = find(abs(BaseTime1-curTime)<1e-12);
            endBase = find(abs(BaseTime2-curTime)<1e-12);
            BaseMFI = mean(curMFI(startBase:endBase));
            plot([BaseTime1,BaseTime2],[BaseMFI,BaseMFI],'-k')
            [PeakTime1,~] = selectPoint(curTime,curMFI);
            [PeakTime2,~] = selectPoint(curTime,curMFI);
            startPeak = find(abs(PeakTime1-curTime)<1e-12);
            endPeak = find(abs(PeakTime2-curTime)<1e-12);
            if length(startPeak:endPeak) >= 3
                peakFit = polyfit(curTime(startPeak:endPeak),curMFI(startPeak:endPeak),2);
                dMFI = polyder(peakFit);
                peakTime = roots(dMFI);
                if peakTime > PeakTime1 && peakTime < PeakTime2 && peakFit(1) < 0
                    PeakMFI = polyval(peakFit,peakTime);
                    plot(curTime(startPeak:endPeak),polyval(peakFit,curTime(startPeak:endPeak)),'-k')
                else
                    PeakMFI = mean(curMFI(startPeak:endPeak));
                    plot([PeakTime1,PeakTime2],[PeakMFI,PeakMFI],'-k')
                end
            else
                PeakMFI = mean(curMFI(startPeak:endPeak));
                plot([PeakTime1,PeakTime2],[PeakMFI,PeakMFI],'-k')
            end
            fluxRatio = (PeakMFI-BG) / (BaseMFI-BG);
            close(gcf)
        else
            BaseMFI = alignData{j+1,12};
            PeakMFI = alignData{j+1,13};
            fluxRatio = alignData{j+1,9};
        end
        peakTimes = cell2mat(alignData(j+1,4:8)) - zeroTime;
        peakTimes = peakTimes(peakTimes > 0);
        peaks = [peakTimes', interp1(curTime,curMFI,peakTimes')];
        
        MFIStored{j} = [curTime, curMFI];
        startInterpTime = max([ceil(curTime(1)/2)*2,-50]);
        endInterpTime = min([floor(curTime(end)/2)*2,200]);
        curTimeInterp = startInterpTime:2:endInterpTime;
        MFIInterp = interp1q(curTime,curMFI,curTimeInterp');
        for k = 1:length(MFIInterp)
            [~,curIdx] = min(abs(timeInterp-curTimeInterp(k)));
            sortedMFI{curIdx} = [sortedMFI{curIdx},MFIInterp(k)];
        end
        % arrange all the data in a sensible way
        if ~isfield(MFIDataSummary,'loc')
            MFIDataSummary(1).date = curDate;
            MFIDataSummary(1).cellNum = curCellNum;
            MFIDataSummary(1).absTime = timeData;
            MFIDataSummary(1).relTime = timeData - zeroTime;
            MFIDataSummary(1).loc = MFIData(curCellNum).loc;
            MFIDataSummary(1).MFI = MFIData(curCellNum).MFI;
%             MFIDataSummary(1).Ca = 3.5*140*(MFIDataSummary(1).MFI-.6)./(2.1-MFIDataSummary(1).MFI);
            MFIDataSummary(1).Ca = 3.892*140*(MFIDataSummary(1).MFI-.552)./(2.1484-MFIDataSummary(1).MFI);
            MFIDataSummary(1).peaks = peaks;
            MFIDataSummary(1).baseline = BaseMFI;
            MFIDataSummary(1).initPeak = PeakMFI;
        else
            Lcur = length(MFIDataSummary)+1;
            MFIDataSummary(Lcur).date = curDate;
            MFIDataSummary(Lcur).cellNum = curCellNum;
            MFIDataSummary(Lcur).absTime = timeData;
            MFIDataSummary(Lcur).relTime = timeData - zeroTime;
            MFIDataSummary(Lcur).loc = MFIData(curCellNum).loc;
            MFIDataSummary(Lcur).MFI = MFIData(curCellNum).MFI;
%             MFIDataSummary(Lcur).Ca = 3.5*140*(MFIDataSummary(Lcur).MFI-.6)./(2.1-MFIDataSummary(Lcur).MFI);
            MFIDataSummary(Lcur).Ca = 3.892*140*(MFIDataSummary(Lcur).MFI-.552)./(2.1484-MFIDataSummary(Lcur).MFI);
            MFIDataSummary(Lcur).peaks = peaks;
            MFIDataSummary(Lcur).baseline = BaseMFI;
            MFIDataSummary(Lcur).initPeak = PeakMFI;
        end
    end
    loadMoreQuest = questdlg('Load more?','Load more?','Yes','No, done','Yes');
    if strcmp(loadMoreQuest,'No, done')
        loadMore = false;
    end
end

%%
timeInterp = -50:2:200;
sortedMFI = cell(1,length(timeInterp));
% updateRatios = questdlg('Update baselines and peaks for new ratios?','Update ratios?','Yes','No','No');
figure
hold on
for i = 1:length(MFIDataSummary)
    curMFI = MFIDataSummary(i).MFI;
    selPts = ~isnan(curMFI);
    curMFI = curMFI(selPts);
    curDate = MFIDataSummary(i).date;
    if ~ischar(curDate)
        curDate = datestr(MFIDataSummary(i).date);
    end
%     if contains(curDate,'2021')
        curMFI = 3.892*140*(curMFI-.552)./(2.1484-curMFI); % convert to Ca
%     elseif contains(curDate,'2022')
%         curMFI = 10*140*(curMFI-.5)./(5-curMFI); % convert to Ca
%         curMFI = 6*140*(curMFI-.5)./(3-curMFI); % convert to Ca
%     else
%         error('What year is it anyway?')
%     end
%     if strcmp(updateRatios,'Yes')
%         cBase = MFIDataSummary(i).baseline;
%         ratioBase = (cBase*2.1 + 140*2.1)/(cBase + 140*2.1/.6);
%         MFIDataSummary(i).baseline = 3.892*140*(ratioBase-.552)./(2.1484-ratioBase);
%         cPeak = MFIDataSummary(i).initPeak;
%         ratioPeak = (cPeak*2.1 + 140*2.1)/(cPeak + 140*2.1/.6);
%         MFIDataSummary(i).initPeak = 3.892*140*(ratioPeak-.552)./(2.1484-ratioPeak);
%     end
    curTime = MFIDataSummary(i).relTime;
    curTime = curTime(selPts);
    [curMFI_noOL,logicOL] = rmoutliers(curMFI,'movmed',9);
    curTime = curTime(~logicOL);
    curMFI = curMFI_noOL;
%     curMFI = interp1(curTime(~logicOL),curMFI_noOL,curTime);
    plot(curTime,curMFI)
    
    startInterpTime = max([ceil(curTime(1)/2)*2,-50]);
    endInterpTime = min([floor(curTime(end)/2)*2,200]);
    curTimeInterp = startInterpTime:2:endInterpTime;
    MFIInterp = interp1(curTime,curMFI,curTimeInterp');
    for k = 1:length(MFIInterp)
        [~,curIdx] = min(abs(timeInterp-curTimeInterp(k)));
        sortedMFI{curIdx} = [sortedMFI{curIdx},MFIInterp(k)];
    end
end
%%
meanMFI = zeros(1,length(sortedMFI));
stdMFI = zeros(1,length(sortedMFI));
semMFI = zeros(1,length(sortedMFI));
for i = 1:length(sortedMFI)
    curSet = sortedMFI{i};
    curSet = curSet(curSet > 0 & curSet < 1000);
    meanMFI(i) = mean(curSet);
    stdMFI(i) = std(curSet);
    semMFI(i) = std(curSet) / sqrt(length(curSet));
end
% plot(timeInterp,meanMFI,'LineWidth',2)
errorbar(timeInterp,meanMFI,semMFI,'LineWidth',1)

%% pull out overall metrics from MFIDataSummary struct
figure
hold on
auc = zeros(length(MFIDataSummary),1);
numPks = zeros(1,length(MFIDataSummary));
aucTimeRange = 100;
for i = 1:length(MFIDataSummary)
    curRange = ~isnan(MFIDataSummary(i).MFI);
    curTime = MFIDataSummary(i).relTime(curRange);
    curCa = MFIDataSummary(i).Ca(curRange);
    [curCa_noOL,logicOL] = rmoutliers(curCa,'movmed',9);
    curCa = interp1(curTime(~logicOL),curCa_noOL,curTime);
    plot(curTime,curCa)
    if curTime(end) < aucTimeRange
        auc(i) = nan;
    else
        startInt = find(curTime > 0, 1, 'first');
        endInt = find(curTime < aucTimeRange, 1, 'last');
        intTime = [0; curTime(startInt:endInt); aucTimeRange];
        intCa = [interp1(curTime,curCa,0); curCa(startInt:endInt); interp1(curTime,curCa,aucTimeRange)];
        auc(i) = trapz(intTime, intCa);
    end
%     numPks(i) = size(MFIDataSummary(i).peaks,1);
    numPks(i) = length(findpeaks(curCa,'MinPeakHeight',150,'MinPeakProminence',20));
end

%% solve for optimal ratios
high0 = 1.8;
low0 = .6;
lowCa = 140*(high0/low0)*(.806-low0)/(high0-.806);
highCa = 140*(high0/low0)*(1.17-low0)/(high0-1.17);
fun = @(x) rootHL(x,lowCa,highCa);
x0 = [2.1,.6];
x = fsolve(fun,x0)


%% solve ratios optimization
low1 = 0.4:.01:0.8;
high1 = 1.8:.01:6;
[low1,high1] = meshgrid(low1,high1);
low1 = low1(:);
high1 = high1(:);
low2 = zeros(size(low1));
high2 = zeros(size(high1));
for i = 1:length(low1)
    lowCa = 140*(high1(i)/low1(i))*(.806-low1(i))/(high1(i)-.806);
    highCa = 140*(high1(i)/low1(i))*(1.17-low1(i))/(high1(i)-1.17);
    fun = @(x) rootHL(x,lowCa,highCa);
    x0 = [high1(i),low1(i)];
    x = fsolve(fun,x0);
    high2(i) = x(1);
    low2(i) = x(2);
end

%% summarize lyse data
timeSample = -50:20:800;
highLyseAligned = cell(size(timeSample));
zeroLyseAligned = cell(size(timeSample));
for i = 1:length(highLyse)
    interpRatios = interp1(highLyse{i}(:,1),highLyse{i}(:,2),timeSample);
    for j = 1:length(timeSample)
        if ~isnan(interpRatios(j))
            highLyseAligned{j} = [highLyseAligned{j},interpRatios(j)];
        end
    end
end

for i = 1:length(zeroLyse)
    interpRatios = interp1(zeroLyse{i}(:,1),zeroLyse{i}(:,2),timeSample);
    for j = 1:length(timeSample)
        if ~isnan(interpRatios(j))
            zeroLyseAligned{j} = [zeroLyseAligned{j},interpRatios(j)];
        end
    end
end

meanHighLyse = zeros(size(timeSample));
stdHighLyse = zeros(size(timeSample));
meanZeroLyse = zeros(size(timeSample));
stdZeroLyse = zeros(size(timeSample));
for i = 1:length(timeSample)
    meanHighLyse(i) = mean(highLyseAligned{i});
    stdHighLyse(i) = std(highLyseAligned{i});
    meanZeroLyse(i) = mean(zeroLyseAligned{i});
    stdZeroLyse(i) = std(zeroLyseAligned{i});
end
figure
hold on
errorbar(timeSample,smooth(meanHighLyse),stdHighLyse,'LineWidth',1)
errorbar(timeSample,meanZeroLyse,stdZeroLyse,'LineWidth',1)

%% measure "final ratios"
figure
highLyseFinal = zeros(size(highLyse));
for i = 1:length(highLyse)
    curTime = highLyse{i}(:,1);
    curRatios = highLyse{i}(:,2);
    plot(curTime,curRatios)
    hold on
    uiwait(msgbox('Choose start and end of the final ratios for this test'))
    [FlatTime1,~] = selectPoint(curTime,curRatios);
    [FlatTime2,~] = selectPoint(curTime,curRatios);
    startFlat = find(abs(FlatTime1-curTime)<1e-12);
    endFlat = find(abs(FlatTime2-curTime)<1e-12);
    highLyseFinal(i) = mean(curRatios(startFlat:endFlat));
    hold off
end

zeroLyseFinal = zeros(size(zeroLyse));
for i = 1:length(zeroLyse)
    curTime = zeroLyse{i}(:,1);
    curRatios = zeroLyse{i}(:,2);
    plot(curTime,curRatios)
    hold on
    uiwait(msgbox('Choose start and end of the final ratios for this test'))
    [FlatTime1,~] = selectPoint(curTime,curRatios);
    [FlatTime2,~] = selectPoint(curTime,curRatios);
    startFlat = find(abs(FlatTime1-curTime)<1e-12);
    endFlat = find(abs(FlatTime2-curTime)<1e-12);
    zeroLyseFinal(i) = mean(curRatios(startFlat:endFlat));
    hold off
end


function F = rootHL(x,lowCa,highCa)
    % x(1) = High ratio
    % x(2) = Low ratio
    F(1) = 140*(x(1)/x(2)) * (0.74-x(2))/(x(1)-0.74) - lowCa;
    F(2) = 140*(x(1)/x(2)) * (1.357-x(2))/(x(1)-1.357) - highCa;
end
