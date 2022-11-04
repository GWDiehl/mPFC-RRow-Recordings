function TSD_Split = splitRRowTSDs_V2(TSD,SplitVar,LapStart,LapEnd,varargin)

% Takes one or more input TSDs and restricts the data based on any of a
% series of relevent RRow contingencies (Lap type, Site, etc.)
%
% Input
% TSD - Single TSD or cell array of TSDs that you wish to restrict based on
%       behavior (TS also ok).
% SplitVar - Cell array of 1 - N different restriction criteria for the
%       TSD. Multiple entries apply AND logic
% LapStart/End - Start or End time of each RRow lap

% Optional Inputs
% vt - Tracking data structure for use with some splitting variables
% LapData - Used for certain splitting variables fullTS - List of full time
% values to be utalized in the case that
%       splitting/restriction is based on TSD data
% maskData - Logical indicating if the restriction should be via masking
%       the data (nan outside of window of interest) or a full restriction
%       where unwanted data are removed. TS objects always remove unwanted
%       data.

% GWD March 2020

%%
vt = [];
fullTS = [];
LapData = [];
maskData = 1;

process_varargin(varargin);

%%

% If only one TSD and it came in not in a cell
if ~iscell(TSD)
    TSD_Split = {TSD};
    breakCell = 1;
else
    TSD_Split = TSD;
    breakCell = 0;
end
if ~iscell(SplitVar)
    SplitVar = {SplitVar};
end

nTSD = length(TSD_Split);

nZones = size(LapStart,2);
WaitZones = 1:nZones;
OfferZones = WaitZones+nZones;

if ~isempty(LapData)
    [LapStart, LapEnd, LapSite, LapAdvance, LapPreStart,...
        OZEnter,OZExit, WZEnter,WZExit, LZEnter,LZExit, AZEnter,AZExit]...
        = identifyRRowLapTimes(LapData,OfferZones,WaitZones,'squeezeD1',0);
end




% Loop through each of the splitting variables sequentially
for iD = 1:length(SplitVar)
    
    limiterType = 'Time';
    
    caseStart = LapStart;
    caseEnd = LapEnd;
    LapsToUse = true(size(caseStart));
    
    if contains(SplitVar{iD},'Delay')
        splittingData = str2num(strrep(SplitVar{iD},'Delay',''));    
        SplitVar{iD} = 'Delay';
    elseif contains(SplitVar{iD},'Value')
        splittingData = str2num(strrep(SplitVar{iD},'Value',''));
        SplitVar{iD} = 'Value';
    elseif contains(SplitVar{iD},'logIdPhi')
        splittingData = str2num(strrep(SplitVar{iD},'logIdPhi',''));
        SplitVar{iD} = 'logIdPhi';
    elseif contains(SplitVar{iD},'Random')
        splittingData = str2num(strrep(SplitVar{iD},'Random',''));    
        SplitVar{iD} = 'Random';
    elseif contains(SplitVar{iD},'SessTime')
        splittingData = str2num(strrep(SplitVar{iD},'SessTime',''));    
        SplitVar{iD} = 'SessTime';  
    elseif contains(SplitVar{iD},'Site')
        splittingData = str2num(strrep(SplitVar{iD},'Site',''));
        SplitVar{iD} = 'Site';
    elseif contains(SplitVar{iD},'WZStart')
        splittingData = str2num(strrep(SplitVar{iD},'WZStart',''));
        SplitVar{iD} = 'WZStart';
    elseif contains(SplitVar{iD},'WZEnd')
        splittingData = str2num(strrep(SplitVar{iD},'WZEnd',''));
        SplitVar{iD} = 'WZEnd';
    elseif contains(SplitVar{iD},'WZFrac')
        splittingData = str2num(strrep(SplitVar{iD},'WZFrac',''));
        SplitVar{iD} = 'WZFrac';
    elseif contains(SplitVar{iD},'LZStart')
        splittingData = str2num(strrep(SplitVar{iD},'LZStart',''));
        SplitVar{iD} = 'LZStart';
    elseif contains(SplitVar{iD},'LZEnd')
        splittingData = str2num(strrep(SplitVar{iD},'LZEnd',''));
        SplitVar{iD} = 'LZEnd';
    elseif contains(SplitVar{iD},'LZFrac')
        splittingData = str2num(strrep(SplitVar{iD},'LZFrac',''));
        SplitVar{iD} = 'LZFrac';
    end
    
    % Grab what is needed for this particular iteration
    switch SplitVar{iD}
        case {'All' 'None' ''}
            % Nothing needs doing
            continue
            
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
            
            % Feeder Site
        case {'Site1' 'Site2' 'Site3' 'Site4'}
            LapsToUse = LapData.ZoneNum(:,WaitZones)==str2double(SplitVar{iD}(5));
            
            % Section of the track
        case 'OfferZone'
            caseStart = OZEnter;
            caseEnd = OZExit;
        case 'WaitZone'
            caseStart = WZEnter;
            % Earn reward or exit the site, whichever came first
            caseEnd = WZExit;
        case 'LingerZone'
            % Will only be valid on earn laps
            caseStart = LZEnter;
            caseEnd = LZExit;
        case 'AdvanceZone' 
            caseStart = AZEnter;
            
        case 'OfferZone+Pre'
            caseStart = LapPreStart;
            caseEnd = OZExit;

        case {'Reward' 'Rail'}
            restrict_TSD =  RRow_TSD_Shell('RadiusRelative',LapStart,LapEnd,...
                'vt',vt,'LapData',LapData,'fullTS',fullTS);
            if strcmp(SplitVar{iD},'Reward')
                window = [.4 1.5];
            elseif strcmp(SplitVar{iD},'Rail')
                window = [-.3, .25];
            end
            limiterType = 'Data';
            
        case {'Delay' 'Value' 'logIdPhi' 'Site'}
            switch SplitVar{iD}
                case 'Delay'
                    currData = LapData.ZoneDelay(:,OfferZones);
                case 'Value'
                    currData = LapData.Value_H(:,OfferZones);
                case 'logIdPhi'
                    currData = log(LapData.IdPhi(:,OfferZones));
                case 'Site'
                    currData = LapData.ZoneNum(:,WaitZones);
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
            
        case 'SessTime'
            if length(splittingData)==2
                sessStart = min(LapData.EnteringZoneTime(:));
                LapsToUse = caseStart >= splittingData(1)+sessStart & caseStart < splittingData(2)+sessStart;
            else
                error('You want a part of the session but did not identify when')
            end
            
        case {'WZStart' 'WZEnd' 'WZFrac' 'LZStart' 'LZEnd'  'LZFrac'}
            switch SplitVar{iD}
                case {'WZStart' 'WZEnd' 'WZFrac'}
                    caseStart = WZEnter;
                    % Earn reward or exit the site, whichever came first
                    caseEnd = WZExit;
                case {'LZStart' 'LZEnd'  'LZFrac'}
                    caseStart = LZEnter;
                    caseEnd = LZExit;
            end     
            
            startBound = caseStart;
            endBound = caseEnd;
            
            if length(splittingData)==2
                if contains(SplitVar{iD},'Start')
                    caseStart = caseStart + splittingData(1);
                    caseEnd = caseStart + splittingData(2);
                elseif contains(SplitVar{iD},'End')
                    caseStart = caseEnd + splittingData(1);
                    caseEnd = caseEnd + splittingData(2);
                elseif contains(SplitVar{iD},'Frac')
                    eventDur = caseEnd - caseStart;
                    caseStart = caseStart + eventDur*splittingData(1);
                    caseEnd = caseStart + eventDur*splittingData(2);
                end
                % Rebound data based on the original extent to ensure we
                % remain within the period of interest (and dont cross into
                % another unwanted event times). 
                caseStart = max(startBound,caseStart);
                caseEnd = min(endBound,caseEnd);
                
                % NOTE: ONE COULD EXPAND OUTSIDE OF THE GIVEN EVENT WINDOW
                % BUT THERE BECOMES RISK OF IDENTIFYING TIMES ACROSS
                % EVENTS. THIS PROBABLY SHOULD BE DONE IN A DIFFERENT WAY
                % TREATING EACH EVENT COMPLETLY SEPERATE (NOT A SINGLUAR
                % TSD)
                
            else
                error('You want a part of the events but did not identify when')
            end
            
            
        otherwise
            error('Case not yet defined')
    end
    
    if strcmp(limiterType,'Time')
        
        SubStart = caseStart(LapsToUse);
        SubEnd = caseEnd(LapsToUse);
        
        % Both times in the pair are valid
        validTimes = ~isnan(SubStart + SubEnd);
        SubStart = SubStart(validTimes);
        SubEnd = SubEnd(validTimes);
    end
    
    % Now that we have the times/data for restriction, go and actually cut
    % out the desired part of our input TSD/TS
    for iTSD = 1:nTSD
        switch limiterType
            case 'Data'
                if iscell(TSD_Split{iTSD})
                    TSD_Split{iTSD} = cellfun(@(x) restrictdata_CrossTSD(x,restrict_TSD,window(1),window(2),'maskData',maskData),TSD_Split{iTSD},'Uniformoutput',0);
                else
                    TSD_Split{iTSD} = restrictdata_CrossTSD(TSD_Split{iTSD},restrict_TSD,window(1),window(2),'maskData',maskData);
                end
            case 'Time'
                
                if iscell(TSD_Split{iTSD}(1))
                    for iE = 1:numel(TSD_Split{iTSD})
                        if strcmp(class(TSD_Split{iTSD}{iE}),'ts') || ~maskData % ts only restrict
                            TSD_Split{iTSD}{iE} = TSD_Split{iTSD}{iE}.restrict(SubStart,SubEnd);
                        else % tsd/ctsd that should be masked
                            TSD_Split{iTSD}{iE} = TSD_Split{iTSD}{iE}.mask(SubStart,SubEnd,0);
                        end
                    end
                else
                    if strcmp(class(TSD_Split{iTSD}),'ts') || ~maskData % ts only restrict
                        TSD_Split{iTSD} = TSD_Split{iTSD}.restrict(SubStart,SubEnd);
                    else % tsd/ctsd that should be masked
                        TSD_Split{iTSD} = TSD_Split{iTSD}.mask(SubStart,SubEnd,0);
                    end
                    
                end
        end
    end
end

if breakCell
    TSD_Split = TSD_Split{:};
end