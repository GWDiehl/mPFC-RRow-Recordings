function fh = PlotCellExampleResponse(exampleCells,events,rawData,dataType,varargin)


[~,minProportion,~,ResultID,OccID,VarID] = ResultsFieldNames(dataType);

plotType = 'wSEM';
normMethod = 'None';

useShaddedError = 1;
yRange = [];
xRange = [];
cRange = [];

timeLabel = [];
timeWindowLimit = [];

[rootDir, DataDefName,BehavDir] = ID_GWD_Dirs;

process_varargin(varargin);


%% General PETH responses across the full population

% No error values calculated for MI
if ismember(dataType,{'MI','MI_TW'})
    plotType = 'AvgRate'; 
    yRoot = 'MutualInfo';
else
    yRoot = 'AvgRate';
end

exampleSSNs = cellfun(@(x) x(1:15),exampleCells,'UniformOutput',0);
nCells = length(exampleCells);
nEvents = length(events);
legendEntry = cell(nEvents,1);
dataIdx = nan(2,nCells);

if nEvents <=7
    % Get the default plotting colors
    h = figure;
    plotColors = get(gca,'colororder'); close(h)
else
    plotColors = jet(nEvents);
end

% Get any behavioral Task Data that may be needed
[FullLapData, FullSessionData, FullIdPhi] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{exampleSSNs});
waitZones = find(ismember(FullLapData.ZoneID(1,1,:),'WaitZone'));
offerZones = find(ismember(FullLapData.ZoneID(1,1,:),'OfferZone'));

behavSSN = cellfun(@(x,y) [x,'-',y],FullSessionData.RatID,FullSessionData.SessDate,'UniformOutput',0);
nLaps = nan(nCells,nEvents);
for iC = 1:nCells
    % Identifiy how may events occured for each event/cell
    behavIdx = ismember(behavSSN,exampleSSNs{iC});
    for iE = 1:nEvents
        switch events{iE}
            case {'OZEntry' 'OZExit'}
                nLaps(iC,iE) = sum(~isnan(reshape(FullLapData.AcceptOffer(behavIdx,:,offerZones),[],1)));
                legendEntry{iE} = 'All';
            case {'OZEntry_Accept' 'OZExit_Accept' 'WZEntry'}
                nLaps(iC,iE) = nansum(reshape(FullLapData.AcceptOffer(behavIdx,:,offerZones),[],1));
                legendEntry{iE} = 'Accept';
            case {'OZEntry_Skip' 'OZExit_Skip'}
                nLaps(iC,iE) = nansum(reshape(FullLapData.SkipOffer(behavIdx,:,offerZones),[],1));
                legendEntry{iE} = 'Skip';
            case {'OZEntry_Site1' 'OZExit_Site1'}
                nLaps(iC,iE) = sum(~isnan(FullLapData.AcceptOffer(behavIdx,:,offerZones(1))));
                legendEntry{iE} = 'Site1';
            case {'OZEntry_Site2' 'OZExit_Site2'}
                nLaps(iC,iE) = sum(~isnan(FullLapData.AcceptOffer(behavIdx,:,offerZones(2))));
                legendEntry{iE} = 'Site2';
            case {'OZEntry_Site3' 'OZExit_Site3'}
                nLaps(iC,iE) = sum(~isnan(FullLapData.AcceptOffer(behavIdx,:,offerZones(3))));
                legendEntry{iE} = 'Site3';
            case {'OZEntry_Site4' 'OZExit_Site4'}
                nLaps(iC,iE) = sum(~isnan(FullLapData.AcceptOffer(behavIdx,:,offerZones(4))));
                legendEntry{iE} = 'Site4';
        end
    end
    
    % Find which particular cell is required
    dataIdx(1,iC) = find(ismember(rawData.Params.SSN,exampleSSNs{iC}));
    dataIdx(2,iC) = find(ismember(rawData.Params.CellID{dataIdx(1,iC)},exampleCells{iC}));
end

if isfield(rawData.Params,'TimingType')
    warpedTime = strcmp(rawData.Params.TimingType,'WarpedTime');
else
    warpedTime = isfield(rawData.Params,'nTimeSteps');
end

% Get the time bins for plotting
if warpedTime
    if isfield(rawData.Params,'nTimeSteps')
        tempSteps = rawData.Params.nTimeSteps;
    else
        tempSteps = rawData.Params.TimeWindow;
    end
    timeEdges = 0:sum(tempSteps);
    transitionBins = cumsum(tempSteps(1:end-1));
    
    if isempty(timeLabel)
        timeLabel = 'Relative Time';
    end
else
    timeEdges = rawData.Params.TimeWindow(1):rawData.Params.dt:rawData.Params.TimeWindow(2);
    transitionBins = 0;
    if isempty(timeLabel)
        timeLabel = 'Time (sec)';
    end
end
timeCenters = computeBinCenters(timeEdges);

if isempty(timeWindowLimit)
    curTime = true(size(timeCenters));
elseif warpedTime
    curTime = timeWindowLimit;
elseif ~warpedTime
    curTime = timeCenters >= timeWindowLimit(1) & timeCenters <= timeWindowLimit(2);
end
timeCenters = timeCenters(curTime);

[subA,subB] = computeSubplotSize(nCells);

fh = figure;
for iC = 1:nCells
    subplot(subA,subB,iC)
    hold on
    for iE = 1:nEvents
        if ismember(dataType,{'TC'})
            paramsIdx = ismember(rawData.Params.(FieldID),events{iE});
            varNames = length(rawData.Params.TuningCurves{paramsIdx});
            dataCenters = cellfun(@(x) computeBinCenters(x),rawData.Params.BinEdges{paramsIdx},'UniformOutput',0);
            
            if length(varNames) == 1
                timeCenters = dataCenters{1};
            elseif length(varNames) == 2 && ~strcmp(plotType,'Image')
                error('You have to use image plotting for 2D TCs')
            else
                error('You cant plot a TC w/ more than 2 dimensions')                
            end
        end
        
        rateData = rawData.Results.(events{iE}).(ResultID){dataIdx(1,iC)}(dataIdx(2,iC),:);
        % Only pull the error data if it is needed (in which case it must
        % exist)
        if ismember(plotType,{'wSEM' 'wStd'})
            rateStd = rawData.Results.(events{iE}).(VarID){dataIdx(1,iC)}(dataIdx(2,iC),:);
        else
            rateStd = nan(size(rateData));
        end
        
        if ~isempty(OccID)
            avgOcc = rawData.Results.(events{iE}).(OccID){dataIdx(1,iC)}';
        else
            avgOcc = ones(size(rateData));
        end
        rateData(avgOcc<minProportion) = nan;
               
        switch normMethod
            case 'None'
                % Nothing
            case 'ZScore'
                rateData = (rateData - nanmean(rateData(:)))./nanstd(rateData(:));
            case 'MeanSubtract'
                rateData = rateData - nanmean(rateData(:));
            case 'ShuffMeanSubtract'
                shuffMean = rawData.Results.Shuffle.(events{iE}).Avg{dataIdx(1,iC)}(dataIdx(2,iC),:);
                rateData = rateData - shuffMean;
            case 'ShuffZScore'
                shuffMean = rawData.Results.Shuffle.(events{iE}).Avg(dataIdx(2,iC),:);
                shuffStd = rawData.Results.Shuffle.(events{iE}).Std(dataIdx(2,iC),:);
                rateData = (rateData - shuffMean)./shuffStd;
        end
        
        % Trim the data time as needed
        if ~isempty(timeWindowLimit)
            rateData = rateData(:,timeWindowLimit);
            rateStd = rateStd(:,timeWindowLimit);
            avgOcc = avgOcc(:,timeWindowLimit);
        end
        
        switch plotType
            case 'AvgRate'
                plot(timeCenters,rateData)
            case 'wSEM'                
                lapEvents = avgOcc*nLaps(iC,iE);
                if useShaddedError
                    ShadedErrorbar(timeCenters,rateData,rateStd./sqrt(lapEvents),'color',plotColors(iE,:),'lineWidth',1,'marker','none');
                else
                    errorbar(timeCenters,rateData,rateStd./sqrt(lapEvents))
                end
                
            case 'wStd'
                if useShaddedError
                    ShadedErrorbar(timeCenters,rateData,rateStd,'color',plotColors(iE,:),'lineWidth',1,'marker','none');
                else
                    errorbar(timeCenters,rateData,rateStd)
                end
            case 'Image'
                s = imagesc(dataCenters{1},dataCenters{2},rateData');
                set(s,'AlphaData',~isnan(rateData'));
        end
    end   
    
    if ~isempty(yRange)
        if iscell(yRange)
            ylim(yRange{iC})
        else
            ylim(yRange)
        end
    end
    if ~isempty(xRange)
        xlim(xRange)
    end
    if ~isempty(cRange) && strcmp(plotType,'Image')
        caxis(cRange);
    end
    if ~strcmp(plotType,'Image') && ismember(normMethod,{'ZScore' 'MeanSubtract' 'ShuffMeanSubtract' 'ShuffZScore'})
        plotHorizLine(0,{'--k'})
    end
    
    if ismember(dataType,{'TC'})
        xlabel(varNames{1})
        if length(varNames) ==2
            ylabel(varNames{2})
        end
    else
        if strcmp(normMethod,'None')
            ylabel(yRoot)
        else
            ylabel({yRoot;normMethod})
        end
        
            plotVertLine(transitionBins,{'--k'})
    
        xlabel(timeLabel)
        if iC == 1
            legend(legendEntry)
        end
    end
    
    title(exampleCells{iC})    
end


