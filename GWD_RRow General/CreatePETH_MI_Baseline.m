function [Results,Params] = CreatePETH_MI_Baseline(varargin)


% Prelim Directoires & DataDef
[rootDir, DataDefName] = ID_GWD_Dirs;

dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

%
timingType = 'RawTime'; % RawTime   WarpedTime
analysisType = 'RawFR'; % RawFR     MutualInfo

RRowEvent = {'OZEntry' 'OZExit' 'WZEntry' 'WZExit'};
EventLimits = {'' '' '' ''};
timeWindow = [-2 2];
dt = .1;
eventTimeBounds = [-1 1];
boundOverride = [];


MI_Params = struct;
MI_Params.Var = 'Accept';
MI_Params.BehavBins = [nan nan 5]; % Start/End/nBins
MI_Params.SpkBins = [1 99 5]; % StartPrctile/EndPrctile/nBins
MI_Params.UseLapData = 1;


nShuffles = 0; % How many shuffles of phys data do we want
shuffMethod = 'TimeShift';
chunkStep = 5;

process_varargin(varargin);

%%

SSNs = dataDef.SSN(OI);
nSess = length(SSNs);

switch timingType
    case 'RawTime'
        timeEdges = timeWindow(1):dt:timeWindow(2);
        timeCenters = timeEdges(2:end) - diff(timeEdges)/2;
        nTimes = length(timeCenters);
    case 'WarpedTime'
        nTimes = sum(timeWindow);
        timeBinIdx = [0 cumsum(timeWindow)];
        timeCenters = nan(1,nTimes); % Required but not used
    otherwise
        error('Your selected timing field does not exist')
end


nTriggers = length(RRowEvent);
eventID = cell(size(RRowEvent));

if ~iscell(eventTimeBounds)
    eventTimeBounds = repmat({eventTimeBounds},nTriggers,1);
end

%% Initialize Outputs

Results = struct;
for iT = 1:nTriggers
    if isempty(EventLimits{iT})
        eventID(iT) = RRowEvent(iT);
    else
        if iscell(EventLimits{iT})
            tempLimit = cell2mat(cellfun(@(x) ['_',x],EventLimits{iT},'UniformOutput',0));
        else
            tempLimit = ['_',EventLimits{iT}];
        end
        eventID{iT} = [RRowEvent{iT},tempLimit];
    end
    eventID{iT} = replace(eventID{iT},{'.','-',':',','},'_');
    
    Results.(eventID{iT}).AvgResults = cell(nSess,1);
    Results.(eventID{iT}).AvgOcc = cell(nSess,1);
    Results.(eventID{iT}).ResultsVar = cell(nSess,1);
    
    if nShuffles
        Results.Shuffle.(eventID{iT}).Example = cell(nSess,1);
        Results.Shuffle.(eventID{iT}).Avg = cell(nSess,1);
        Results.Shuffle.(eventID{iT}).Std = cell(nSess,1);
    end
end

% Parameters
Params = struct;
Params.SSN = SSNs';
Params.CellID = cell(nSess,1);
Params.CellFR = cell(nSess,1);

Params.TimingType = timingType;
Params.AnalysisType = analysisType;
if strcmp(analysisType,'MutualInfo')
    Params.MI_Params = MI_Params;
end
Params.dt = dt;
Params.TimeWindow = timeWindow;
Params.RRowEvent = RRowEvent;
Params.EventLimits = EventLimits;
Params.FullEventID = eventID;
Params.EventTimeBounds = eventTimeBounds;
Params.BoundsOverride = boundOverride;

Params.nShuffles = nShuffles;
if nShuffles
    Params.ShuffleMethod = shuffMethod;  
end

%% General Behav Info

% Get any behavioral Task Data that may be needed
[FullLapData, FullSessionData] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});

[~,nLaps,nZones] = size(FullLapData.Decision);
waitZones = find(ismember(squeeze(FullLapData.ZoneID(1,1,:)),'WaitZone'));
offerZones = find(ismember(squeeze(FullLapData.ZoneID(1,1,:)),'OfferZone'));

% Grab the relevant Behav Times
[~,~,~,~,~,~,~,~,~,~,~,~,~,eventTime]...
    = identifyRRowLapTimes(FullLapData,offerZones,waitZones);

%%

for iS = 1:nSess
    
    [vt, ~, S, S_FNs] = GWD_LoadRRowSession(SSNs{iS},{'Path' 'Spikes'},'restrictToTaskTime',1);
    
    nCells = length(S_FNs);
    nSpks = cellfun(@(x) length(x.range),S);
    
    sessStart = vt.x.starttime;
    sessEnd = vt.x.endtime;
    
    lapStart = reshape(eventTime.LapStart(iS,:),nLaps,nZones);
    lapEnd = reshape(eventTime.LapEnd(iS,:),nLaps,nZones);
    
    cellFR = nSpks/(sessEnd - sessStart);
    
    Params.CellID{iS} = S_FNs;
    Params.CellFR{iS} = cellFR;
    
    currLapData = subsetStructByIndx(FullLapData,iS,'squeezeD1',1);
    currSessData = subsetStructByIndx(FullSessionData,iS,'squeezeD1',1);
    
    [~,allRRowEvents] = IdentifyAllRRowEvents(currLapData);
        
    % Generate shuffled data if needed
    preparedSpiking = cell(nShuffles+1,1);
    preparedSpiking{1} = S;
    for iP = 1:nShuffles
        preparedSpiking{iP+1} = shuffleSpkData(S,shuffMethod,'tStart',sessStart,'tEnd',sessEnd);
    end    
    
    % nCells x nShuffles+1 length of spike times
    preparedSpiking = cat(1,preparedSpiking{:});
    
    % Get relevant behavioral data   
    if ~isempty(boundOverride)
        boundOverride.Behav = RRow_TSD_Shell(boundOverride.Var,lapStart,lapEnd,'LapData',currLapData);
    end
    if strcmp(analysisType,'MutualInfo')
        [MI_Behav, MI_TSDBins] = RRow_TSD_Shell_V2(MI_Params.Var,lapStart,lapEnd,...
            'vt',vt,'LapData',currLapData,'SessData',currSessData,'fullTS',vt.x.range);
        % Overwrite binning info for MI calculations
        for iX = 1:3
            if isnan(MI_Params.BehavBins(iX))
                MI_Params.BehavBins(iX) = MI_TSDBins(iX);
            end
        end
    end
    
    for iT = 1:nTriggers
        
        %% Get the respective event trigger times
        [eventTimes, ~, limiterMask] = identifyRRowEvents(currLapData,RRowEvent{iT},'limits',EventLimits{iT});
        triggerTimes = selectData(eventTimes,limiterMask);
        % Note that the triggerIdx runs sequentially in the task, Zones x
        % Laps
        triggerIdx = reshape(1:length(limiterMask(:)),nZones,nLaps)';
        
        % Remove any bad trigger times (nans, or wont contain the full time
        % window due to begining/end of session)
        switch timingType
            case 'RawTime'
                badTimes = triggerTimes<=(sessStart-timeEdges(1)) | triggerTimes>=(sessEnd-timeEdges(end));
            case 'WarpedTime'
                badTimes = triggerTimes<=allRRowEvents(abs(min(eventTimeBounds{iT}))) | triggerTimes>=allRRowEvents(end-max(eventTimeBounds{iT}));
        end
        badTimes = badTimes | isnan(triggerTimes);
        
        triggerTimes = triggerTimes(~badTimes);
        triggerIdx = triggerIdx(~badTimes);
        
        % Sort them so things run chronologically not by site
        [triggerTimes,order] = sort(triggerTimes);
        triggerIdx = triggerIdx(order);
        
        nEvents = length(triggerTimes);
        
        if nEvents == 0
            Results.(eventID{iT}).AvgResults{iS} = nan(nCells,nTimes,1);
            Results.(eventID{iT}).AvgOcc{iS} = nan(nTimes,1);
            Results.(eventID{iT}).ResultsVar{iS} = nan(nCells,nTimes,1);            
            if nShuffles
                Results.Shuffle.(eventID{iT}).Example{iS} = nan(nCells,nTimes,1);
                Results.Shuffle.(eventID{iT}).Avg{iS} = nan(nCells,nTimes,1);
                Results.Shuffle.(eventID{iT}).Std{iS} = nan(nCells,nTimes,1);                
            end
            continue
        end
        
        % Identify temporal bins that sit outside of the behavioral window of
        % interest. First for the firing rates and second for which entries
        % will be used in the regression calcualtion
        [outOfRange,~,eventBounds] = IdentifyTimesWithinBounds(triggerTimes,allRRowEvents,timeCenters,eventTimeBounds{iT});
        
        switch timingType
            case 'RawTime'
                sessTimeEdges = timeEdges;
                sessTimeCenters = repmat(timeCenters,nEvents,1);
                
            case 'WarpedTime'
                outOfRange = false(nTimes,nEvents);
                
                if ~isempty(boundOverride)
                    eventBounds = OverrideEventBounds(eventBounds,triggerTimes,boundOverride.Behav,boundOverride.Val,boundOverride.Ref,boundOverride.Dur);
                end
                eventBounds = sort(eventBounds,2);
                sessTimeEdges = nan(nEvents,nTimes+1);
                sessTimeCenters = nan(nEvents,nTimes);
                for iE = 1:nEvents
                    tempSteps = [];
                    for iB = 1:length(timeWindow)
                        tempSteps = cat(2,tempSteps,linspace(eventBounds(iE,iB),eventBounds(iE,iB+1),timeWindow(iB)+1));
                    end
                    sessTimeEdges(iE,:) = unique(tempSteps) - triggerTimes(iE);
                    sessTimeCenters(iE,:) = computeBinCenters(sessTimeEdges(iE,:));
                end                
        end
        
        % Grand output FR (Cells*nShuffles+1 x Time x Events)
        cellFRData = ExtractChunkedFRData(preparedSpiking,triggerTimes,sessTimeEdges,'chunkStep',chunkStep);
        cellFRData(:,outOfRange) = nan;  
                      
        switch analysisType
            case 'RawFR'                
                % (Cells x Time x Events x Shuffles)
                curatedFRData = permute(reshape(cellFRData,nCells,nShuffles+1,nTimes,nEvents),[1 3 4 2]);
                
                % Collect outputs for the session
                Results.(eventID{iT}).AvgResults{iS} = nanmean(curatedFRData(:,:,:,1),3);
                Results.(eventID{iT}).ResultsVar{iS} = nanstd(curatedFRData(:,:,:,1),[],3);
                
                if nShuffles                    
                    Results.Shuffle.(eventID{iT}).Example{iS} = nanmean(curatedFRData(:,:,:,2),3);
                    Results.Shuffle.(eventID{iT}).Avg{iS} = nanmean(nanmean(curatedFRData(:,:,:,2:end),3),4);
                    Results.Shuffle.(eventID{iT}).Std{iS} = nanstd(nanmean(curatedFRData(:,:,:,2:end),3),[],4);
                end
                
            case 'MutualInfo'
                fullBehav_MI = nan(nTimes,nEvents);
                for iE = 1:nEvents
                    if MI_Params.UseLapData
                        fullBehav_MI(:,iE) = MI_Behav.data(triggerTimes(iE));
                    else
                        fullBehav_MI(:,iE) = MI_Behav.data(triggerTimes(iE)+sessTimeCenters(iE,:));
                    end
                end
                fullBehav_MI(outOfRange) = nan;
                
                % Discretize the behavior data for MI calculations
                fullBehav_MI = discretize(fullBehav_MI,linspace(MI_Params.BehavBins(1),MI_Params.BehavBins(2),MI_Params.BehavBins(3)+1));
                
                % Verify things are sequential integers so MI does not
                validBins = ~isnan(fullBehav_MI);
                [~,~,temp] = unique(fullBehav_MI(validBins));
                fullBehav_MI(validBins) = temp;
                
                infoResults = nan(nCells*(nShuffles+1),nTimes);
                % Run calculations
                for iC = 1:nCells*(nShuffles+1)
                    % Discretize the spiking data
                    tempCell = reshape(cellFRData(iC,:),nTimes,nEvents);
                    tempCell = discretize(tempCell,linspace(prctile(tempCell(:),MI_Params.SpkBins(1)),prctile(tempCell(:),MI_Params.SpkBins(2)),MI_Params.SpkBins(3)+1));
                    
                    % Verify things are sequential integers so MI does not
                    validBins = ~isnan(tempCell);
                    [~,~,temp] = unique(tempCell(validBins));
                    tempCell(validBins) = temp;
                    
                    for iX = 1:nTimes
                        inputData = cat(2,reshape(fullBehav_MI(iX,:),[],1),reshape(tempCell(iX,:),[],1));
                        
                        % You only have one behavioral element; MI is not
                        % valid
                        if max(inputData(:,1)) == 1
                            continue
                        end
                        infoResults(iC,iX) = calculateInfoMetric(inputData,'PairedMI',0); % Dont check 'unique'
                    end
                end
                
                % Collect outputs for the session
                % (Cells x Time x Shuffles)
                curatedInfoData = permute(reshape(infoResults,nCells,nShuffles+1,nTimes),[1 3 2]);
                Results.(eventID{iT}).AvgResults{iS} = curatedInfoData(:,:,1);
                
                if nShuffles
                    Results.Shuffle.(eventID{iT}).Example{iS} = curatedInfoData(:,:,2);
                    Results.Shuffle.(eventID{iT}).Avg{iS} = nanmean(curatedInfoData(:,:,2:end),3);
                    Results.Shuffle.(eventID{iT}).Std{iS} = nanstd(curatedInfoData(:,:,2:end),[],3);
                end
        end
        Results.(eventID{iT}).AvgOcc{iS} = nanmean(~outOfRange,2);        

    end
    
    fprintf('Done with %s: Session %d of %d \n',SSNs{iS},iS,nSess)
end