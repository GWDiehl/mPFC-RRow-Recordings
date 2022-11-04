function [fh_Sig] = PlotSingleCellSignificance(sigTime,Params,dataCenters,anatGroup,varargin)

selectedTime = [8:33];
GWDColor = loadGWD_RRowColormap;
xRange = [];
plotTitle = '';

process_varargin(varargin);

if isempty(xRange)
    xRange = [selectedTime(1)-1.5 selectedTime(end)+.5];
end

if isfield(Params,'TimingType')
    warpedTime = strcmp(Params.TimingType,'WarpedTime');
else
    warpedTime = isfield(Params,'nTimeSteps');
end

if warpedTime
    if isfield(Params,'nTimeSteps')
        tempSteps = Params.nTimeSteps;
    else
        tempSteps = Params.TimeWindow;
    end
    transitionBins = cumsum(tempSteps(1:end-1));
else
    transitionBins = 0;
end

nAnat = length(anatGroup);
offSetSize = .5/nAnat/2;

nPlots = size(sigTime,3);
[subA,subB] = computeSubplotSize(nPlots);


fh_Sig = figure;
for iP = 1:nPlots
    subplot(subA,subB,iP)
    hold on
    for iA = 1:nAnat
        plot(dataCenters,sigTime(iA,:,iP)+(iA-.5-nAnat/2)*offSetSize,'o','color',GWDColor.(anatGroup{iA}))
    end
    
    xlim(xRange)
    ylim([-.2 1.2])
    plotVertLine(transitionBins,{'--k'})
    
    xlabel('Relative Time');
    ylabel('Significant Data')
    
    if iP == 1
        legend(anatGroup,'location','best')
    end
    
    if iscell(plotTitle)
        title(plotTitle{iP})
    else
        title(plotTitle)
    end
end

