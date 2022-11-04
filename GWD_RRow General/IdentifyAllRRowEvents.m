function [allEvents,sortedEvents] = IdentifyAllRRowEvents(behavLapData)


switch size(behavLapData.ZoneID,3)
    case 1        
        waitZones = find(ismember(behavLapData.ZoneID(1,:),'WaitZone'));
        offerZones = find(ismember(behavLapData.ZoneID(1,:),'OfferZone'));
        fullData = 0;
    otherwise
        waitZones = find(ismember(behavLapData.ZoneID(1,1,:),'WaitZone'));
        offerZones = find(ismember(behavLapData.ZoneID(1,1,:),'OfferZone'));
        fullData = 1;
end


[~, ~, ~, ~, ~,OZEnter,~,WZEnter,~,LZEnter,~,AZEnter,~]...
    = identifyRRowLapTimes(behavLapData,offerZones,waitZones);


% Gather up all RRow zone transition times, remove nans, and sort by
% time for use in limiting to within a given task state window.
allEvents = cat(1+fullData,OZEnter,WZEnter,LZEnter,AZEnter);

if size(allEvents,3)>1
    sortedEvents = [];
else
    sortedEvents = sort(reshape(allEvents(~isnan(allEvents)),[],1));
end