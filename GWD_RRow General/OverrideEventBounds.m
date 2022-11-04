
function adaptedBounds = OverrideEventBounds(eventBounds,eventTimes,overrideData,overrideVal,overrideRef,overrideDur);

% If we want to override some of the values in our event boundaries
% according to specific behavioral components around that particular event
% time. Check if each event qualifies and then override the event boundary
% time according to the selected duration around the selected boundary
% event.

% GWD May 2022


adaptedBounds = eventBounds;
[nEvents,nBounds] = size(eventBounds);

for iE = 1:nEvents
    if overrideData.data(eventTimes(iE)) >= overrideVal(1) && overrideData.data(eventTimes(iE)) <= overrideVal(2)
        for iB = 1:nBounds
            if ~isnan(overrideDur(iB))
                adaptedBounds(iE,iB) = eventBounds(iE,overrideRef(iB)) + overrideDur(iB);
            end
        end
    end
end