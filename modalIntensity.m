function modeVal = modalIntensity(im)
    pixelVals = double(im(:));
    pixelVals = pixelVals(pixelVals > 0);
    xInt = std(pixelVals)/100;
%     [xVec,density,xMax] = kernelDensity(xInt,pixelVals,'parabolic');
    xVec = min(pixelVals):xInt:max(pixelVals);
    [density,xVec] = ksdensity(pixelVals,xVec);
%     histObj = histogram(pixelVals,500);
%     xVals = histObj.BinEdges;
%     count = histObj.Values;
%     xGap = histObj.BinWidth;
    [count,xVals] = histcounts(pixelVals,100);
    xGap = mean(diff(xVals));
%     if xGap < 2
%         xGap = 2;
%         histObj = histogram(pixelVals,'BinWidth',xGap);
%     end
    areaUnderHist = sum(count) * xGap;
    densityPlot = density * areaUnderHist;
%     hold on
%     plot(xVec,densityPlot,'LineWidth',2)
    [~,maxIdx] = max(density);
    % use golden algorithm to find max
    gr = (1+sqrt(5))/2;
    a = xVec(maxIdx-1);
    b = xVec(maxIdx+1);
    c = b - (b-a)/gr;
    d = a + (b-a)/gr;
    tol = 0.01;
    while abs(c-d) > tol
        [densTests,~] = ksdensity(pixelVals,[c,d]);
        if densTests(1) > densTests(2)
            b = d;
        else
            a = c;
        end
        c = b - (b-a)/gr;
        d = a + (b-a)/gr;
    end
    xPeak = (a+b)/2;
%     yL = get(gca,'ylim');
%     plot([xPeak, xPeak], yL, 'Color', 'k', 'LineWidth', 2)
%     xlim([xPeak-xInt*200, xPeak+xInt*200])
    modeVal = xPeak;
end