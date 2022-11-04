function [allEvents,sortedEvents,eventIdx,sortedIdx] = IdentifySelectedRRowEvents(behavLapData,selectedEvents)


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


[~,~,~,~,~,~,~,~,~,~,~,~,~,FullEventTimes] = ...
    identifyRRowLapTimes(behavLapData,offerZones,waitZones);

nEvents = length(selectedEvents);


% Gather up all RRow zone transition times, remove nans, and sort by
% time for use in limiting to within a given task state window.
allEvents = FullEventTimes.(selectedEvents{1});
% Identify which of the set of event selections each event time corresponds
% to
eventIdx = repmat(1,size(allEvents));
for iX = 2:nEvents
    allEvents = cat(1+fullData,allEvents,FullEventTimes.(selectedEvents{iX}));
    eventIdx = cat(1+fullData,eventIdx,repmat(iX,size(FullEventTimes.(selectedEvents{iX}))));
end

if size(allEvents,3)>1
    sortedEvents = [];
    sortedIdx = [];
else
    [sortedEvents,order] = sort(reshape(allEvents(~isnan(allEvents)),[],1));
    temp = reshape(eventIdx(~isnan(allEvents)),[],1);
    sortedIdx = temp(order);
end