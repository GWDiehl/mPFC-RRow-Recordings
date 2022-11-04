function [Results, Params] = RRow_InfoTheory_TimeWarpped(varargin)


% Computes RRow MI calculations analysis using a time warped depiction of
% the tast event. Things are treated not as raw time but instead as
% proportion of the time in a particular zone.
%
%
%
% OUPTUS
% Results - A Structure containing fields with each PETH
%   calculated. Each has a cell array of 1 x NSess with each cell
%   containing nCells x nBins mtx of Avg Firing rate and Rate Variation
%   (std).
% Params - General set of parameters of what was calculated, which
%   sessions/cells are included, and what settings were used.
%
% Optional Inputs to override defaults:
% RRowEvent - Cell array of which Trigger Events you would like to align to.
% EventLimits - Matching cell array of any limiters that you would like to
%       restrict the events by (eg. Site Number, Choice).
% NOTE -- RRowEvent and EventLimits will be combined into EventID which
%       will define the results structure field. Find these combos in
%       Params.
% limitEventBounds - From the RRowEvent of interest how many events do we
%       want to be looking at within the sequence. Note that all of the
%       event times/task periods will be concatinated and sorted in
%       chronoloical order.
% nTimeSteps - How many time steps do we want for each task event window
%       included in the full set of events.
% overrideTSD - A TSD of potential by lap/time to trigger an override of
%       the full event bounds.
% overrideValue - Window of values for which the TSD will trigger an
%       override
% overrideDuration - If we are overriding how much time should be used in
%       the task window (nan will bypass overriding for particular event
%       windows)
% overrideRef - Which event entry in the original bounds are we refrencing
%       against w/ overrideDruation to produce the override.
%
% nShuffles - How many shuffles of the data do you want to compute along
%   with the real data. 0 corresponds to no shuffles, only real.
% shuffMethod - What method of shuffling do you want to use?


% GWD May 2021

%%


% Prelim Directoires & DataDef
[rootDir, DataDefName] = ID_GWD_Dirs;

dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

%Analysis Method
Method = 'PairedMI';

% Behavior Variables
behavTC = 'Accept';
useFullTS = 0;
nBins = 5;
TCBins_Override = [];

RRowEvent = {'OZEntry'};
EventLimits = {''};
subsetParam = [];

limitEventBounds = [-1 0 1 2];
nTimeSteps = [20 40 20];

overrideTSD = 'Accept';
overrideValue = [.5 1.5];
overrideDuration = [nan nan nan 1];
overrideRef = [nan nan nan 3];

nShuffles = 0; % How many shuffles of phys data do we want
shuffMethod = 'TimeShift';
commonTShift = 0;

chunkStep = 5;

process_varargin(varargin);

%%

if ~ismember(Method,{'PairedMI' 'Entropy'})
    error('Unknown Information Theory Method')
end

nTimes = sum(nTimeSteps);

subDataDef = subsetFromDataDef(dataDef,OI);

SSNs = subDataDef.SSN;
nSess = length(SSNs);

% Get any behavioral Task Data that may be needed
[FullLapData, FullSessionData] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});

nZones = size(FullLapData.Decision,3);
WaitZones = 1:nZones;
OfferZones = WaitZones+nZones;

nTriggers = length(RRowEvent);
eventID = cell(size(RRowEvent));

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
    
    Results.(eventID{iT}).Info = cell(nSess,1);
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

Params.RRowEvent = RRowEvent;
Params.EventLimits = EventLimits;
Params.FullEventID = eventID;
Params.EventSubsampling = subsetParam;

Params.TuningCurve = behavTC;
if ~isempty(TCBins_Override)
    % Check and log any TCBin Overrides
    assert(length(TCBins_Override) ==2,'Your Bin Override input is not a 2 element vector')
    Params.TCBins_Override = TCBins_Override;
end

Params.nTimeSteps = nTimeSteps;
Params.LimitEventBounds = limitEventBounds;

Params.OverrideTSD = overrideTSD;
Params.OverrideValue = overrideValue;
Params.OverrideDuration = overrideDuration;

Params.nShuffles = nShuffles;
if nShuffles
    Params.ShuffleMethod = shuffMethod;
    Params.CommonShift = commonTShift;
end

for iS = 1:nSess
    
    [vt, ~, S, S_FNs] = GWD_LoadRRowSession(SSNs{iS},{'Path' 'Spikes'},'restrictToTaskTime',1);
    
    if useFullTS
        fullTS = vt.x.range;
    else
        fullTS = [];
    end
        
    nCells = length(S_FNs);
    nSpks = cellfun(@(x) length(x.range),S);
    
    sessStart = vt.x.starttime;
    sessEnd = vt.x.endtime;
    
    cellFR = nSpks/(sessEnd - sessStart);
    
    Params.CellID{iS} = S_FNs;
    Params.CellFR{iS} = cellFR;
            
    % Generate shuffled data if needed
    preparedSpiking = cell(nShuffles+1,1);
    preparedSpiking{1} = S;
    for iP = 1:nShuffles
        preparedSpiking{iP+1} = shuffleSpkData(S,shuffMethod,'tStart',sessStart,'tEnd',sessEnd,'commonShiftAcrossCells',commonTShift);
    end
    
     %% Get the behavoral data
    currLapData = subsetStructByIndx(FullLapData,iS,'squeezeD1',1);
    currSessData = subsetStructByIndx(FullSessionData,iS,'squeezeD1',1);
    
    [LapStart, LapEnd, LapSite, LapAdvance, LapPreStart,...
        OZEnter,OZExit, WZEnter,WZExit, LZEnter,LZExit, AZEnter,AZExit]...
        = identifyRRowLapTimes(currLapData,OfferZones,WaitZones);    
    
    % Gather up all RRow zone transition times, remove nans, and sort by
    % time for use in limiting to within a given task state window.
    allRRowEvents = cat(1,OZEnter,WZEnter,LZEnter,AZEnter);
    allRRowEvents = allRRowEvents(~isnan(allRRowEvents));
    allRRowEvents = sort(allRRowEvents(:));
    
    if ~isempty(overrideTSD)
        TSD_Override = RRow_TSD_Shell_V2(overrideTSD,LapStart,LapEnd,'LapData',currLapData);
    end
    
    [BehavTSD, TCBins] = RRow_TSD_Shell_V2(behavTC,LapStart,LapEnd,...
        'bins',nBins,'vt',vt,'LapData',currLapData,'SessData',currSessData,'fullTS',fullTS);
    
    % Limit the range of your TC Variable as anything outside of the bin
    % edges will discretizes as a NaN
    if ~isempty(TCBins_Override)
        TCBins(1:2) = TCBins_Override;
    end
    
    Params.BehaviorBins = TCBins(3);
    binEdges = linspace(TCBins(1),TCBins(2),TCBins(3)+1);
    BehavTSD.D = discretize(BehavTSD.D,binEdges);
    
    for iZ = 1:nTriggers
        
        %% Get the respective event trigger times
        [eventTimes, ~, limiterMask] = identifyRRowEvents(currLapData,RRowEvent{iZ},'limits',EventLimits{iZ});
        triggerTimes = selectData(eventTimes,limiterMask);
        
        % Subset the trigger times according to additional criteria of
        % desired
        if iscell(subsetParam)
            triggerTimes = MatchRRowEventCount(triggerTimes,subsetParam{iZ},currLapData);
        end
        
        % Remove any nans
        triggerTimes = triggerTimes(~isnan(triggerTimes));
        
        % Enseure that we only treat trigger times that have enough outside
        % events
        if min(limitEventBounds) < 0
            triggerTimes(triggerTimes<=allRRowEvents(abs(min(limitEventBounds)))) = [];
        end
        if max(limitEventBounds) > 0
            triggerTimes(triggerTimes>=allRRowEvents(end-max(limitEventBounds))) = [];
        end      
        
        % Sort them so things run chronologically not by site
        triggerTimes = sort(triggerTimes);
        
        nEvents = length(triggerTimes);
                        
        if nEvents <2 % Must have at least two events to calculate Info Metrics
            % Collect outputs
            if strcmp(Method,'Entropy')
                Results.(eventID{iZ}).Info{iS} = nan(nCells+1,nTimes);
            else
                Results.(eventID{iZ}).Info{iS} = nan(nCells,nTimes);
            end
            
            if nShuffles
                Results.Shuffle.(eventID{iZ}).Example{iS} = nan(nCells,nTimes);
                Results.Shuffle.(eventID{iZ}).Avg{iS} = nan(nCells,nTimes);
                Results.Shuffle.(eventID{iZ}).Std{iS} = nan(nCells,nTimes);
            end
            continue
        end
        
        surroundBounds = nan(nEvents,length(limitEventBounds));
        
        for iE = 1:nEvents
            for iB = 1:length(limitEventBounds)
                if limitEventBounds(iB) < 0
                    tempIdx = FastFind(allRRowEvents - triggerTimes(iE) + .01 < 0,abs(limitEventBounds(iB)),'last');
                    if isempty(tempIdx)
                        surroundBounds(iE,iB) = triggerTimes(iE) - 1;
                    else
                        surroundBounds(iE,iB) = allRRowEvents(tempIdx(1));
                    end
                elseif limitEventBounds(iB) > 0
                    tempIdx = FastFind(allRRowEvents - triggerTimes(iE) - .01 > 0,abs(limitEventBounds(iB)),'first');
                    if isempty(tempIdx)
                        surroundBounds(iE,iB) = triggerTimes(iE) + 1;
                    else
                        surroundBounds(iE,iB) = allRRowEvents(tempIdx(end));
                    end
                elseif limitEventBounds(iB) == 0
                    surroundBounds(iE,iB) = triggerTimes(iE);
                end
            end
            % If we meet a particular override condition for this event
            % dont use the full event window but compute and internally
            % referenced override time within.
            if ~isempty(overrideTSD)
                if TSD_Override.data(triggerTimes(iE)) >= overrideValue(1) && TSD_Override.data(triggerTimes(iE)) <= overrideValue(2)
                    origBounds = surroundBounds(iE,:);
                    for iB = 1:length(limitEventBounds)
                        if ~isnan(overrideDuration(iB))
                            surroundBounds(iE,iB) = origBounds(overrideRef(iB)) + overrideDuration(iB);
                        end
                    end
                end
            end
        end
        
       % Make sure that for each event the times run in sequential order
        fullEventTimes = sort(surroundBounds,2);
        fullTimeSteps = nan(nEvents,sum(nTimeSteps)+1);
        
        behavData = nan(nTimes,nEvents);
        for iE = 1:nEvents
            tempSteps = [];
            for iB = 1:length(nTimeSteps)
                tempSteps = cat(2,tempSteps,linspace(fullEventTimes(iE,iB),fullEventTimes(iE,iB+1),nTimeSteps(iB)+1));
            end
            fullTimeSteps(iE,:) = unique(tempSteps);
            
            % Are we pulling a lap based TSD (entries ~= nLaps)?
            % If so generalize off of the trigger time as that is the
            % lap data of interest.
            if length(BehavTSD.range) < sum(~isnan(LapStart(:)))*2.5
                behavData(:,iE) = BehavTSD.data(triggerTimes(iE));
            else
                behavData(:,iE) = BehavTSD.data(computeBinsCenters(fullTimeSteps));
            end
        end
        
        % Collect All the spiking (Real + Shuffles)
        completeSpks = cat(1,preparedSpiking{:});
        % Grand output FR (Cells x Shuff x Time x Events)
        cellData = nan(nCells,nShuffles+1,nTimes,nEvents);
                
        if chunkStep
            chunkedEvents = nan(chunkStep,ceil(nEvents/chunkStep));
            chunkedEvents(1:nEvents) = 1:nEvents;
            
            for iE = 1:size(chunkedEvents,1)
                tempIdx = chunkedEvents(iE,:);
                tempIdx(isnan(tempIdx)) = [];
                if isempty(tempIdx)
                    continue
                end
                QTimes = fullTimeSteps(tempIdx,:)';
                if ~issorted(QTimes(:))
                    error('You have overlap in your timing, use a larger chunk size');
                end
                QTemp = MakeQfromS_NonUniformTime(completeSpks, ts(QTimes(:)),'tStart',sessStart,'tEnd',sessEnd);
                QTemp = recenterQMtx(QTemp);
                for iF = 1:length(tempIdx)
                    triggIdx = tempIdx(iF);
                    currEventTimes = computeBinCenters(fullTimeSteps(triggIdx,:));
                    spikeCounts = full(QTemp.data(currEventTimes))';
                    spikeRate = spikeCounts./diff(fullTimeSteps(triggIdx,:));
                    % Reshape the Shuff back to D2
                    cellData(:,:,:,triggIdx) = reshape(spikeRate,nCells,nShuffles+1,sum(nTimeSteps));
                end
            end
        else            
            % Deal with each event entriley independently
            for iE = 1:nEvents
                fullTimeSteps = [];
                for iB = 1:length(nTimeSteps)
                    fullTimeSteps = cat(2,fullTimeSteps,linspace(fullEventTimes(iE,iB),fullEventTimes(iE,iB+1),nTimeSteps(iB)+1));
                end
                
                fullTimeSteps = unique(fullTimeSteps);
                Q = MakeQfromS_NonUniformTime(completeSpks,ts(fullTimeSteps),'tStart',sessStart,'tEnd',sessEnd);
                spikeCounts = full(Q.D(1:end-1,:))';
                spikeRate = spikeCounts./diff(fullTimeSteps);
                
                % Reshape the Shuff back to D2
                cellData(:,:,:,iE) = reshape(spikeRate,nCells,nShuffles+1,sum(nTimeSteps));
            end
        end
        
        % Convert FRs into binned data (use only current real data FRs to
        % find bounds)
        binnedCellData = nan(size(cellData));
        for iC = 1:nCells
            
            FR_Bounds = prctile(cellData(iC,1,:),[1.5 98.5]);
            tempCellData = cellData(iC,:);
            tempBinned = discretize(tempCellData,linspace(FR_Bounds(1),FR_Bounds(2),nBins+1));
            
            % Fill in any FRs that fell outside of the percentiles
            % (currently just ignore/NaN these outliers)
% %             tempBinned(tempCellData<FR_Bounds(1)) = 1;
% %             tempBinned(tempCellData>FR_Bounds(2)) = bins;            
            binnedCellData(iC,:) = tempBinned;            
        end
        
        InfoResults = nan(nCells,nTimes,nShuffles+1);
        % Run calculations
        for iC = 1:nCells  % Parfor gets an indexing crash for no good reason. Does not seem to work so STOP TRYING.
            for iX = 1:nShuffles+1
                tempCellData = binnedCellData(iC,iX,:,:);
                for iT = 1:nTimes
                    switch Method
                        case 'PairedMI'
                            inputData = cat(2,behavData(iT,:)',squeeze(tempCellData(:,:,iT,:)));
                        case 'Entropy'
                            inputData = squeeze(tempCellData(:,:,iT,:));
                    end
                    InfoResults(iC,iT,iX) = calculateInfoMetric(inputData,Method);
                end
            end
        end
                        
        if strcmp(Method,'Entropy') % Get the entropy for the behavior
            behavEntropy = nan(1,nTimes,nShuffles+1);
            for iT = 1:nTimes
                inputData = behavData(iT,:)';
                behavEntropy(1,iT,:) = calculateInfoMetric(inputData,Method);
            end
            InfoResults = cat(1,InfoResults,behavEntropy);
        end
        
        % Collect outputs
        Results.(eventID{iZ}).Info{iS} = InfoResults(:,:,1);
        
        if nShuffles
            Results.Shuffle.(eventID{iZ}).Example{iS} = InfoResults(:,:,2);
            Results.Shuffle.(eventID{iZ}).Avg{iS} = nanmean(InfoResults(:,:,2:end),3);
            Results.Shuffle.(eventID{iZ}).Std{iS} = nanstd(InfoResults(:,:,2:end),[],3);
        end
        
    end
    
    fprintf('Done with analysis for %s, %d of %d \n',SSNs{iS},iS,nSess)
    
end




