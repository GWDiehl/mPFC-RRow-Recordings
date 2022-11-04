function [Results,Params] = ExtractAvgRateBetweenEvents(varargin)


% Prelim Directoires & DataDef
[rootDir, DataDefName] = ID_GWD_Dirs;

dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');


eventStart = {'LapStart' 'OZEntry' 'WZEntry' 'LZEntry'};
eventEnd = {'LapEnd' 'OZExit' 'WZExit' 'LZExit'};
EventLimits = {'' '' '' ''};

nShuffles = 5; % How many shuffles of phys data do we want
shuffMethod = 'TimeShift';
chunkStep = 5;

process_varargin(varargin);

%%

SSNs = dataDef.SSN(OI);
nSess = length(SSNs);

nTriggers = length(eventStart);
eventID = cell(size(eventStart));


%% Initialize Outputs

Results = struct;
for iT = 1:nTriggers
    if isempty(EventLimits{iT})
        eventID{iT} = [eventStart{iT},'_',eventEnd{iT}];
    else
        if iscell(EventLimits{iT})
            tempLimit = cell2mat(cellfun(@(x) ['_',x],EventLimits{iT},'UniformOutput',0));
        else
            tempLimit = ['_',EventLimits{iT}];
        end
        eventID{iT} = [eventStart{iT},'_',eventEnd{iT},tempLimit];
    end
    % Add any additional suffix that has been passed in
    eventID{iT} = replace(eventID{iT},{'.','-',':',','},'_');
    
    Results.(eventID{iT}).AvgResults = cell(nSess,1);
    Results.(eventID{iT}).LapIdx = cell(nSess,1);
    
    if nShuffles
        Results.Shuffle.(eventID{iT}).Example = cell(nSess,1);
        Results.Shuffle.(eventID{iT}).Avg = cell(nSess,1);
        Results.Shuffle.(eventID{iT}).Std = cell(nSess,1);
    end
end

assert(isunique(eventID),['Your event entries are not all unique. ',...
    'The structure/output will not keep them distinct so add in suffix fields'])

% Parameters
Params = struct;
Params.SSN = SSNs';
Params.CellID = cell(nSess,1);
Params.CellFR = cell(nSess,1);

Params.EventStart = eventStart;
Params.EventEnd = eventEnd;
Params.FullEventID = eventID;

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
    
    for iT = 1:nTriggers
        
        %% Get the respective event trigger times
        [startTimes, ~, mask_Start] = identifyRRowEvents(currLapData,eventStart{iT},'limits',EventLimits{iT});
        [endTimes, ~, mask_End] = identifyRRowEvents(currLapData,eventEnd{iT},'limits',EventLimits{iT});
        fullMask = mask_Start & mask_End;
        
        startTimes = selectData(startTimes,fullMask);
        endTimes = selectData(endTimes,fullMask);
        
        % Note that the triggerIdx runs sequentially in the task, Zones x
        % Laps
        triggerIdx = reshape(1:length(fullMask(:)),nZones,nLaps)';
        
        % Remove any bad trigger times (nans, or wont contain the full time
        % window due to begining/end of session)
        badTimes = isnan(startTimes) | isnan(endTimes);
        
        startTimes = startTimes(~badTimes);
        endTimes = endTimes(~badTimes);
        triggerIdx = triggerIdx(~badTimes);
        
        % Sort them so things run chronologically not by site
        [startTimes,order] = sort(startTimes);
        endTimes = endTimes(order);
        triggerIdx = triggerIdx(order);
        
        nEvents = length(startTimes);
        
        % Combine start and end times to a single bounds variable, and take
        % the middle (avg) as the "Trigger Time" for purposes of utalizing
        % existing code.
        eventBounds = cat(2,startTimes,endTimes);
        triggerTimes = nanmean(eventBounds,2);
        eventBounds = eventBounds-triggerTimes;
        
        %%% This approch of making a Q matrix of non-uniform bins and
        %%% chunking is about 3-5x faster than restricting each
        %%% cell around each event (coded above). Results are identical
        %%% (TESTED and VERIFIED)
        % Grand output FR (Cells*nShuffles+1 x Time x Events)        
        cellFRData = ExtractChunkedFRData(preparedSpiking,triggerTimes,eventBounds,'chunkStep',chunkStep);
        
        % Cells*nShuffles+1 x Events
        cellFRData = permute(cellFRData,[1 3 2]);
        
        curatedFRData = permute(reshape(cellFRData,nCells,nShuffles+1,nEvents),[1 3 2]);
        
        % Collect outputs for the session
        Results.(eventID{iT}).AvgResults{iS} = curatedFRData(:,:,1);
        Results.(eventID{iT}).LapIdx{iS} = triggerIdx;
        
        if nShuffles
            Results.Shuffle.(eventID{iT}).Example{iS} = curatedFRData(:,:,2);
            Results.Shuffle.(eventID{iT}).Avg{iS} = nanmean(curatedFRData(:,:,2:end),3);
            Results.Shuffle.(eventID{iT}).Std{iS} = nanstd(curatedFRData(:,:,2:end),[],3);
        end
        
    end
    fprintf('Done with %s: Session %d of %d \n',SSNs{iS},iS,nSess)
end
