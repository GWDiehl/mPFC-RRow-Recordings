function [fh, outputResults, dataCenters, Params] = PlotFullSingleCellResults(dataFN,dataType,selectedEvent,selectedGroup,plotTitle,anatGroup,varargin)

% Plot multiple groups of mPFC single cell responding during RRow on the
% same plot. Plots anatGroups all on the same plot and selected Event/Group
% on different subplots of the same figure.

% GWDiehl Oct 2021


GWDColor = loadGWD_RRowColormap;
selectedTime = [8:33];
windowBuffer = [0 0];
yRange = [];
plotType = 'ShaddedErrorBar';
timeLabel = '';
normMethod = 'ShuffMeanSubtract';

anatGrpFN = 'E:\DATA-Geoff\GDiehl Data\Analysis\CellAnatGroupings';
groupingLogic = 'AND';

process_varargin(varargin);

%%

switch dataType
    case 'PETH'
        dataLabel = 'Normalized FR';
    case 'MI'
        dataLabel = 'Normalzied Mutual Info';
    otherwise
        dataLabel = '';
end

nAnat = length(anatGroup);
nPlots = length(selectedEvent);
[subA,subB] = computeSubplotSize(nPlots);
outputResults = cell(nPlots,nAnat);

%
fh = figure;
for iP = 1:nPlots
    
    [selectedData,cellIDs,Params] = ExtractPFCAnalysisData(dataFN,dataType,selectedEvent{iP},'normMethod',normMethod);
    if isfield(Params,'TimingType')
        warpedTime = strcmp(Params.TimingType,'WarpedTime');
    else
        warpedTime = isfield(Params,'nTimeSteps');
    end
        
    % Identify which cells for this plotting set and identify the
    % anatomical location of each
    selectedCells = SelectPFCSubset(cellIDs,selectedGroup{iP},'logicType',groupingLogic);
    
    if isempty(anatGrpFN)
        anatSelection = groupPFCCellsByDepth(cellIDs,anatGroup);
    else
       anatSelection = AllocatePFCAnatToLocigals(anatGrpFN,cellIDs,anatGroup);
    end
    
    
    % Get the time bins for plotting
    if warpedTime
        if isfield(Params,'nTimeSteps')
            tempSteps = Params.nTimeSteps;
        else
            tempSteps = Params.TimeWindow;
        end
        binEdges = 0:sum(tempSteps);
        transitionBins = cumsum(tempSteps(1:end-1));
            
        if isempty(timeLabel)
            timeLabel = 'Relative Time';
        end
    else
        binEdges = Params.TimeWindow(1):Params.dt:Params.TimeWindow(2);
        transitionBins = 0;
        if isempty(timeLabel)
            timeLabel = 'Time (sec)';
        end
    end
    dataCenters = computeBinCenters(binEdges);
    
    if isempty(selectedTime)
        curTime = true(size(dataCenters));
    elseif warpedTime
        curTime = selectedTime;
    elseif ~warpedTime
        curTime = dataCenters >= selectedTime(1) & dataCenters <= selectedTime(2);
    end
    dataCenters = dataCenters(curTime);
    xRange = [dataCenters(1) + windowBuffer(1),dataCenters(end) + windowBuffer(2)];
    
    subplot(subA,subB,iP); hold on
    
    for iA = 1:length(anatGroup)                
        if sum(selectedCells & anatSelection(:,iA)) == 0
            continue
        end
        outputResults{iP,iA} = selectedData(selectedCells & anatSelection(:,iA),curTime);
        PlotSingleCellAvg(outputResults{iP,iA},dataCenters,'plotType',plotType,'plotColor',GWDColor.(anatGroup{iA}))
    end
    
    if ~isempty(yRange)
        ylim(yRange)
    end
    if ~isempty(xRange)
        xlim(xRange)
    else
        xlim([binEdges(1),binEdges(end)])
    end
    
    plotVertLine(transitionBins,{'--k'})
    plotHorizLine(0,{'--k'})
    
    xlabel(timeLabel);
    ylabel(dataLabel)
    if iP == 1 && ~strcmp(plotType,'ShaddedErrorBar')
        legend(anatGroup,'location','best')
    end    
    title(plotTitle{iP})
end

