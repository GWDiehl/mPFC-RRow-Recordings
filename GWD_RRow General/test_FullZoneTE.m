
function [Results, Params] = test_FullZoneTE(varargin)


% Prelim Directoires & DataDef
[rootDir, DataDefName] = ID_GWD_Dirs;

dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

RRowEvent = {'OfferZone' 'OfferZone' 'OfferZone' 'WaitZone' 'LingerZone'};
EventLimits = {'' 'Skip' 'Accept' '' '' ''};

dt = .1;

nShuffles = 5; % How many shuffles of phys data do we want
shuffMethod = 'TimeShift';

fullTS = [];
useFullTS = 0;

validCells = {'All'};

process_varargin(varargin);

chunkOptions = {'OfferZone' 'WaitZone' 'LingerZone' 'AdvanceZone'};

%%

SSNs = dataDef.SSN(OI);
nSess = length(SSNs);

% Get any behavioral Task Data that may be needed
[FullLapData] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});

waitZones = find(ismember(FullLapData.ZoneID(1,1,:),'WaitZone'));
offerZones = find(ismember(FullLapData.ZoneID(1,1,:),'OfferZone'));


%%

nTriggers = length(RRowEvent);
eventID = cell(size(RRowEvent));

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

%% Parameters

Params = struct;
Params.SSN = SSNs';
Params.CellID = cell(nSess,1);
Params.CellFR = cell(nSess,1);
Params.dt = dt;
Params.RRowEvent = RRowEvent;
Params.EventLimits = EventLimits;
Params.FullEventID = eventID;
Params.ValidCells = validCells;

Params.nShuffles = nShuffles;
if nShuffles
    Params.ShuffleMethod = shuffMethod;
end


for iS = 1:nSess
    
    [vt, ~, S, S_FNs] = GWD_LoadRRowSession(SSNs{iS},{'Path' 'Spikes'},'restrictToTaskTime',1);
    
    currLapData = subsetStructByIndx(FullLapData,iS,'squeezeD1',1);
    [LapStart, LapEnd] = identifyRRowLapTimes(currLapData,offerZones,waitZones,'tsdStep',0);
    
    if useFullTS
        fullTS = vt.x.range;
    end
            
    mazeChunkData = RRow_TSD_Shell('MazeChunk',LapStart,LapEnd,'LapData',currLapData,'fullTS',fullTS);
    
    cellsToUse = groupPFCCellsByDepth(S_FNs,validCells);
    cellsToUse = any(cellsToUse,2);
    
    nCells = length(S);
    nSpks = cellfun(@(x) length(x.range),S);
    
    sessStart = vt.x.starttime;
    sessEnd = vt.x.endtime;
    
    cellFR = nSpks/(sessEnd - sessStart);    
    
    Params.CellID{iS} = S_FNs;
    Params.CellFR{iS} = cellFR;
    
    % Generate shuffled data if needed
    preparedSpiking = cell(nShuffles+1,1);
    Q = MakeQfromS(S,dt,'tStart',sessStart,'tEnd',sessEnd);
    preparedSpiking{1} = recenterQMtx(Q);
    for iP = 1:nShuffles
        shuffS = shuffleSpkData(S,shuffMethod,'tStart',sessStart,'tEnd',sessEnd);
        shuffQ = MakeQfromS(shuffS,dt,'tStart',sessStart,'tEnd',sessEnd);
        preparedSpiking{iP+1} = recenterQMtx(shuffQ);
    end
    nSpkBins = length(preparedSpiking{1}.range);
    
    for iT = 1:nTriggers
        % Limit the chunk data according to event limits (Skip/Accept)
        triggerChunkData = splitRRowTSDs_V2(mazeChunkData,EventLimits{iT},LapStart,LapEnd,'LapData',currLapData);
        
        % Identify the specific chunk value for the desired event
        mazeChunkIdx = find(ismember(chunkOptions,RRowEvent{iT}));
        if ismember(RRowEvent{iT},{'FullSess','All'})
            mazeChunkBounds = [0.5 length(chunkOptions)+.5];
        else
            mazeChunkBounds = [mazeChunkIdx-.5,mazeChunkIdx+.5];
        end
        
        % Mask the spiking TSDs according to the behavioral conditions
        currRestrictSpk = preparedSpiking;
        [~,validBehavBins] = restrictdata_CrossTSD(preparedSpiking{1},triggerChunkData,mazeChunkBounds(1),mazeChunkBounds(2));
        validBehavBins = validBehavBins.data;
        
        % Valid times for the TE calculations have current and prior binns
        % included
        validTETimes = validBehavBins(2:end) & validBehavBins(1:end-1);
        
        InfoResults = nan(nCells,nCells,nShuffles+1);
        rawData = full(currRestrictSpk{1}.data);
        % For each cell make the data run sequentially from 1. Do this once now so
        % you can bypass in the TE calculation
        for iC = 1:nCells
            cellData = rawData(:,iC);
            validBins = validBehavBins;
            [~,~,sequentialCell] = unique(cellData(validBins));
            rawData(validBins,iC) = sequentialCell;
        end
        
        % Run calculations
        for iP = 1:nShuffles+1 % First run is the real data and each subsequent is a shuffle
            currShuffData = full(currRestrictSpk{iP}.data);

            % Make the data run sequentially from 1
            for iC = 1:nCells
                cellData = currShuffData(:,iC);
                validBins = validBehavBins;
                [~,~,sequentialCell] = unique(cellData(validBins));
                currShuffData(validBins,iC) = sequentialCell;
            end

            %%% Unsure if this parfor will work, or if it will throw stupid indx errors
            parfor iC = 1:nCells
                if ~cellsToUse(iC)
                    continue
                end
                for iC2 = 1:nCells
                    % Dont do same cell and by def the first time
                    % point wont work (needs lag data)
                    if iC == iC2 || ~cellsToUse(iC2)
                        continue
                    end
                    % [Y-Future, Y-Past, X-Past]
                    % (Note: Only X is shuffled == Random degree of
                    % influence from cell X/iC2 given this specific
                    % configuration of Y responses Past/Future)
                    yFut = rawData(2:end,iC);
                    yPast = rawData(1:end-1,iC);
                    xPast = currShuffData(1:end-1,iC2);
                    
                    inputData = cat(2,yFut(validTETimes),yPast(validTETimes),xPast(validTETimes));
                    
                    % iC == Output; iC2 == input
                    InfoResults(iC,iC2,iP) = calculateInfoMetric(inputData,'TransferEntropy',0);
                end
            end
        end
        
        % Collect outputs
        Results.(eventID{iT}).Info{iS} = InfoResults(:,:,1);
        
        if nShuffles            
            Results.Shuffle.(eventID{iT}).Example{iS} = InfoResults(:,:,2);
            Results.Shuffle.(eventID{iT}).Avg{iS} = nanmean(InfoResults(:,:,2:end),3);
            Results.Shuffle.(eventID{iT}).Std{iS} = nanstd(InfoResults(:,:,2:end),[],3);
        end
    end
    
    fprintf('Done with analysis for %s, %d of %d \n',SSNs{iS},iS,nSess)
    
    
end
