function [splitTimes,dataMask,splitLower,splitUpper] = splitRRowPETH(EventTimes,curSplit,LapData,TriggerEvent,varargin)

% Restricts the input event times for PETH analyses beased on a desired set
% of splits (1 or more). This cut down the orinal EventTimes vector with
% only those that fall within the given split as well as a mask of what
% remained from the original input data and upper and lower time bounds
% around each resultant event time signifying the edges of the restriction
% window (relevent when the PETH window is larger than the restriction
% window; i.e. time could be inculded in the PETH that falls outside of the
% restriction window)
%
% Inputs
% EventTimes - vector of the original PETH event trigger times
% curSplit - cell array (or can be string if only one) of the split
%       variables that are desired. Multiple split varibales will apply
%       though AND logic. Currently code is only built to support time
%       windows, not windows generated from specific sets of TSD data (see
%       splitRRowTSDs_V2 for example w/ TSDs).
% LapData - Full set of behaivoral lap data for the session that provides
%       necessary info about restriction time points.
% TriggerEvent - What event was used for generating the original
%       EventTimes. In some cases this will necessitate variations to the
%       bounds such that the right data are maintained after splitting. See
%       below.
%
% OUTPUTS
% splitTimes - Remaining times after restricitng data
% dataMask - A mask of which of the original data passed through
%       restriction.
% splitLower/Upper - Upper and Lower time bounds around each resultatnt
%       spitTime that indicate the temporal bounds of the closest
%       restriction window. Relevent in cases where PETH window exceeds
%       restriction window.

% GWD May 2020

% Timeframe for the iniital/adaptation segment of a lap (in minutes)
lapInitialSegment = 10; 

process_varargin(varargin);

nZones = size(LapData.Decision,2);
WaitZones = 1:nZones;
OfferZones = WaitZones+nZones;

[LapStart, LapEnd, ~, LapAdvance, LapPreStart] ...
    = identifyRRowLapTimes(LapData,OfferZones,WaitZones);

if ~iscell(curSplit)
    curSplit = {curSplit};
end
SplitDim = length(curSplit);
splittingData = [];

splitTimes = EventTimes(:); % Make sure these are a vector
splitLower = 0;
splitUpper = 0;

for iD = 1:SplitDim
    FullStart = LapStart;
    FullEnd = LapEnd;
    LapsToUse = true(size(FullStart));
    
    % A few events require a bit of special splitting to make sure that
    % they include the relevant time windows as there are inherent clashes
    % between wanted PETH window and splitting window.
    
    % Note: Other variables will completely omit any relevent events (e.g.
    % feeder times on skip laps). The current code will do this as told
    % even though it leaves no data.
    if ismember('OfferZoneEntry',TriggerEvent)
        % Take from the advance zone of the previous lap
        FullStart = LapPreStart;
    elseif ismember('WaitZoneEntry',TriggerEvent)
        % Take from entry into the offer zone
        FullStart = LapData.EnteringZoneTime(:,OfferZones);
        
    elseif ismember('AdvanceZoneExit',TriggerEvent)
        % Take from the exit of the subsequent Offer ZOne
        NextOZExit = advanceRRowData(LapData.ExitZoneTime(:,OfferZones),-1)-1e-10;
        NextOZExit(end) = nan; % Kill the first lap that gets looped to the back
        
        FullEnd = NextOZExit;
    end
    
    if contains(curSplit{iD},'Delay')
        splittingData = str2num(strrep(curSplit{iD},'Delay',''));    
        curSplit{iD} = 'Delay';
    elseif contains(curSplit{iD},'Value')
        splittingData = str2num(strrep(curSplit{iD},'Value',''));
        curSplit{iD} = 'Value';
    elseif contains(curSplit{iD},'logIdPhi')
        splittingData = str2num(strrep(curSplit{iD},'logIdPhi',''));
        curSplit{iD} = 'logIdPhi';
    elseif contains(curSplit{iD},'Random')
        splittingData = str2num(strrep(curSplit{iD},'Random',''));    
        curSplit{iD} = 'Random';
    elseif contains(curSplit{iD},'SessTime')
        splittingData = str2num(strrep(curSplit{iD},'SessTime',''));    
        curSplit{iD} = 'SessTime';
    elseif contains(curSplit{iD},'SessionTime')
        splittingData = str2num(strrep(curSplit{iD},'SessionTime',''));
        curSplit{iD} = 'SessionTime';
    elseif contains(curSplit{iD},'PercentTotalEarns')
         splittingData = str2num(strrep(curSplit{iD},'PercentTotalEarns',''));    
        curSplit{iD} = 'PercentTotalEarns';
    elseif contains(curSplit{iD},'Site')
        splittingData = str2num(strrep(curSplit{iD},'Site',''));
        curSplit{iD} = 'Site';
    elseif contains(curSplit{iD},'Rank')
        splittingData = str2num(strrep(curSplit{iD},'Rank',''));
        curSplit{iD} = 'Rank';
    elseif contains(curSplit{iD},'PrctLingerTime')
        splittingData = str2num(strrep(curSplit{iD},'PrctLingerTime',''));
        curSplit{iD} = 'PrctLingerTime';
    elseif contains(curSplit{iD},'LingerTime')
        splittingData = str2num(strrep(curSplit{iD},'LingerTime',''));
        curSplit{iD} = 'LingerTime';
    end
    
    switch curSplit{iD}
        case 'All'
            % No special changes other than lap start/end
            
            %%% Take only a subset of laps
            
            % Behavioral Choice
        case 'Skip'
            LapsToUse = LapData.SkipOffer(:,OfferZones)==1;
        case 'Accept'
            LapsToUse = LapData.AcceptOffer(:,OfferZones)==1;
        case 'Quit'
            LapsToUse = LapData.QuitOffer(:,WaitZones)==1;
        case 'Earn'
            LapsToUse = LapData.EarnOffer(:,WaitZones)==1;
        case 'NoEarn'
            LapsToUse = LapData.QuitOffer(:,WaitZones)==1 | LapData.SkipOffer(:,OfferZones)==1;
            
            
        case 'SessInitial'
            LapTime = LapData.EnteringZoneTime(:,OfferZones);
            LapTime = LapTime - min(LapTime(:));
            LapsToUse = LapTime < (lapInitialSegment * 60);
        case 'SessPostInitial'
            LapTime = LapData.EnteringZoneTime(:,OfferZones);
            LapTime = LapTime - min(LapTime(:));
            LapsToUse = LapTime > (lapInitialSegment * 60);
            
        case {'Delay' 'Value' 'logIdPhi' 'Site' 'Rank' 'LingerTime' 'PrctLingerTime' 'PercentTotalEarns'}
            switch curSplit{iD}
                case 'Delay'
                    currData = LapData.ZoneDelay(:,OfferZones);
                case 'Value'
                    currData = LapData.Value_H(:,OfferZones);
                case 'logIdPhi'
                    currData = log(LapData.IdPhi(:,OfferZones));
                case 'Site'
                    currData = LapData.ZoneNum(:,WaitZones);
                case 'Rank'
                    currData = LapData.SiteRank(:,WaitZones);
                case {'LingerTime','PrctLingerTime'}
                    currData = ExtractBehavDataFromFull(LapData,[],curSplit{iD});
                case 'PercentTotalEarns'                    
                    currData = ExtractBehavDataFromFull(LapData,[],curSplit{iD});
                    
            end
            switch length(splittingData)
                case 2
                    LapsToUse = currData >= splittingData(1) & currData <= splittingData(2);                    
                case 1
                    LapsToUse = currData == splittingData;
                case 0
                    error('You want a %s but did not specify what kind',curSplit{iD})
                otherwise
                   LapsToUse = ismember(currData,splittingData);
            end
            
            
        case {'SessTime' 'SessionTime'}
            if length(splittingData)==2
                sessStart = min(LapData.EnteringZoneTime(:));
                LapsToUse = EventTimes > splittingData(1)+sessStart & EventTimes < splittingData(2)+sessStart;
            else
                error('You want a part of the session but did not identify when')
            end                    
            
        case 'Random'            
            if length(splittingData) == 1
                if splittingData > 1
                    splittingData = splittingData/100;
                end
                LapsToUse = false(size(FullStart));
                validTimes = ~isnan(EventTimes);
                for iZ = 1:size(validTimes,2)
                    validIdx = find(validTimes(:,iZ));
                    lapsToKeep = randsample(validIdx,round(length(validIdx)*splittingData));
                    LapsToUse(lapsToKeep,iZ) = 1;
                end
            else
                error('You want a random subset but misidentified how much')
            end            
      
            
            %%% Get a new set of Start/End times, these will
            %%% only have the valid times so no need for LapsToUse
            % Section of the track
        case 'OfferZone'
            if ~ismember('OfferZoneEntry',TriggerEvent)
                FullStart = LapData.EnteringZoneTime(:,OfferZones);
            end
            if ~ismember('OfferZoneExit',TriggerEvent)
                FullEnd = LapData.ExitZoneTime(:,OfferZones);
            end
        case 'WaitZone'
            if ~ismember('WaitZoneEntry',TriggerEvent)
                FullStart = LapData.EnteringZoneTime(:,WaitZones);
            end
            if ~ismember('WaitZoneExit',TriggerEvent)
                % Earn reward or exit the site, whichever came first
                FullEnd = min(LapData.ExitZoneTime(:,WaitZones),LapData.FeederTimes(:,WaitZones));
            end
        case 'LingerZone'
            % Will only be valid on earn laps
            FullStart = LapData.FeederTimes(:,WaitZones);
            FullEnd = LapData.ExitZoneTime(:,WaitZones);
        case 'AdvanceZone'
            FullStart = LapAdvance;
        otherwise
            error('Case not yet defined')
    end
    
    SubStart = FullStart(LapsToUse);
    SubEnd = FullEnd(LapsToUse);
    
    
    %%% Trim down times as needed
    
    % Both times in the pair are valid
    validTimes = ~isnan(SubStart + SubEnd);
    SubStart = SubStart(validTimes);
    SubEnd = SubEnd(validTimes);
    
    % Only keep the events within the Splitting Windows
    %%% NOTE: Splitting windows (SubStart/SubEnd) MUST be non-overlapping
    %%% such that an event can only fall within a single splitting window.
    %%% In standrd event bounds this should not be a problem.
    validEvents = iswithin(splitTimes,SubStart,SubEnd);
    nEvents = sum(validEvents);
    splitTimes = splitTimes(validEvents);
    
    previousLower = splitLower;
    previousUpper = splitUpper;
    
    splitLower = nan(nEvents,1);
    splitUpper = nan(nEvents,1);
    
    % Find the bounds that correspond to each Event Time
    for iSE = 1:nEvents
        % Multiple events may share the same bounds so they are
        % repeated for each one and can just use the first
        % instance, however times SHOULD NEVER fall within two
        % distinct sets of bounds.
        prevIdx = FastFind(arrayfun(@(x,y) iswithin(splitTimes(iSE),x,y),previousLower,previousUpper),1);
        curIdx = FastFind(arrayfun(@(x,y) iswithin(splitTimes(iSE),x,y),SubStart,SubEnd),1);
        % Narrow the bounds as needed
        if ~any(prevIdx)
            splitLower(iSE) = SubStart(curIdx);
            splitUpper(iSE) = SubEnd(curIdx);
        else
            splitLower(iSE) = max(previousLower(prevIdx),SubStart(curIdx));
            splitUpper(iSE) = min(previousUpper(prevIdx),SubEnd(curIdx));
        end
    end
    
end

% Which of the original timestamps made it through the splitting
dataMask = ismember(EventTimes,splitTimes);
