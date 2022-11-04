function [eventTimes, limitedTimes, limiterMask] = identifyRRowEvents(lapData,Event,varargin)

limits = [];

process_varargin(varargin);


nZones = size(lapData.Decision,2);
WaitZones = 1:nZones;
OfferZones = WaitZones+nZones;

[LapStart, LapEnd, LapSite, LapAdvance, LapPreStart,...
    OZEnter,OZExit, WZEnter,WZExit, LZEnter,LZExit, AZEnter,AZExit,...
    FullEventTimes, WZExit_Adj, LZExit_Adj, AZEnter_Adj] = identifyRRowLapTimes(lapData,OfferZones,WaitZones);

switch Event
    case 'FullSess'
        eventTimes = nanmean(LapStart(:));
        limits = [];
        
    case fieldnames(FullEventTimes)
        eventTimes = FullEventTimes.(Event);
        
    case 'AllZoneEntry'
        eventTimes = cat(1,OZEnter,WZEnter,LZEnter,AZEnter);
        
    case 'All-OZ'
        eventTimes = cat(1,WZEnter,LZEnter,AZEnter);
    case 'All-WZ'
        eventTimes = cat(1,OZEnter,LZEnter,AZEnter);
    case 'All-LZ'
        eventTimes = cat(1,OZEnter,WZEnter,AZEnter);
    case 'All-AZ'
        eventTimes = cat(1,OZEnter,WZEnter,LZEnter);
        
        % I was not consistant with using Enter vs Entry!!!. Now this has
        % to be here still to catch the differences
    case 'OZEntry'
        eventTimes = OZEnter;
    case 'WZEntry'
        eventTimes = WZEnter;
    case 'LZEntry'
        eventTimes = LZEnter;
    case 'AZEntry'
        eventTimes = AZEnter;
        
    case 'MidOZ'
        eventTimes = (OZEnter+OZExit)/2;
    case 'MidWZ'
        eventTimes = (WZEnter+WZExit)/2;
        
    case {'RandomNLaps' 'RandomAccept' 'RandomSkip'}
        switch Event
            case 'RandomNLaps'
                eventTimes = LapStart';
            case 'RandomAccept'
                eventTimes = selectData(OZExit,lapData.AcceptOffer(:,OfferZones) == 1)';
            case 'RandomSkip'
                eventTimes = selectData(OZExit,lapData.SkipOffer(:,OfferZones) == 1)';
        end
        nEvents = sum(~isnan(eventTimes(:)));
        randTimes = sort(randi(round([min(LapStart(:)),max(LapEnd(:))]),nEvents,1));
        eventTimes(~isnan(eventTimes)) = randTimes;
        eventTimes = eventTimes';
        
    case 'Skip'
        eventTimes = selectData(OZExit,lapData.SkipOffer(:,OfferZones) == 1);
    case 'Accept'
        eventTimes = selectData(OZExit,lapData.AcceptOffer(:,OfferZones) == 1);
    case 'Quit'
        eventTimes = selectData(WZExit,lapData.QuitOffer(:,WaitZones) == 1);
    case 'Earn'
        eventTimes = selectData(WZExit,lapData.EarnOffer(:,WaitZones) == 1);
    case 'PostRewardAdvance'
        eventTimes = selectData(AZEnter,lapData.EarnOffer(:,WaitZones) == 1);
        
    otherwise
        error('Case not yet defined')
end

if ~isempty(limits)
    [limitedTimes, limiterMask]= splitRRowPETH(eventTimes,limits,lapData,Event);
else
    limitedTimes = eventTimes;
    limiterMask = true(size(eventTimes));
end