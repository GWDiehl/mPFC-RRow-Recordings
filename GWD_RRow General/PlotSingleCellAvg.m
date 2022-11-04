function PlotSingleCellAvg(inputData,xBins,varargin)

yRange = [];
xRange = [];
plotType = 'ShaddedErrorBar';
plotColor = [];

centerFun = @nanmean;

process_varargin(varargin);

dataMean = centerFun(inputData,1);
dataStd = nanstd(inputData,[],1);
dataCnt = sum(~isnan(inputData),1);

if isempty(dataMean)
    return
end

switch plotType
    case 'ShaddedErrorBar'
        dataSEM = dataStd./sqrt(dataCnt);
        if isempty(plotColor)
            ShadedErrorbar(xBins,dataMean,dataSEM,'lineWidth',1,'marker','none');
        else
            ShadedErrorbar(xBins,dataMean,dataSEM,'color',plotColor,'lineWidth',1,'marker','none');
        end
    case 'ErrorBar'
        dataSEM = dataStd./sqrt(dataCnt);
        if isempty(plotColor)
            errorbar(xBins,dataMean,dataSEM);
        else
            errorbar(xBins,dataMean,dataSEM,'color',plotColor);
        end
    case 'LinePlot'
        if isempty(plotColor)
            plot(xBins,dataMean)
        else
            plot(xBins,dataMean,'color',plotColor)
        end
end

if ~isempty(yRange)
    ylim(yRange)
end
if ~isempty(xRange)
    xlim(xRange)
end