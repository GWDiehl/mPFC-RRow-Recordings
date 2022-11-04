function [outOfRange, inRange, eventBounds] = IdentifyTimesWithinBounds(selectedEvents,allEvents,timeBins,boundLimits,varargin)

% Identify time points that fall outside of a set of boundary times as
% defined by other behavioral event times. This is directly used in things
% like identifying RRow events that occur only within a particular phase of
% the task and remove any times outside of the window of interest.

% GWDiehl Oct 2021

%%

withinEventRange = 1e-6;
maxBoundRange = inf;

process_varargin(varargin);

%%

if isempty(boundLimits) | boundLimits == 0
    inRange = ~outOfRange;
    return
end

% Identify the number of events and time bins and initialize an output
nEvents = length(selectedEvents);
nTimes = length(timeBins);
outOfRange = false(nTimes,nEvents);

assert(length(boundLimits)>=2,'You need to have at least 2 boundary limits')


allEvents = sort(reshape(allEvents(~isnan(allEvents)),[],1));
eventBounds = nan(nEvents,length(boundLimits));

for iE = 1:nEvents
    tempTimes = selectedEvents(iE) + timeBins;
    
    % Remove the current event from the list of candidate bounaries if it
    % is a member
    alternativeCandidates = allEvents < selectedEvents(iE)-withinEventRange | allEvents > selectedEvents(iE)+withinEventRange;
    candidateBounds = allEvents(alternativeCandidates);
    
    for iB = 1:length(boundLimits)
        
        % Negative bound limits, look for the last events before the event
        % time. Take the first one as that is the requested events away.
        if boundLimits(iB) < 0
            tempIdx = FastFind(candidateBounds - selectedEvents(iE) < 0,abs(boundLimits(iB)),'last');
            if length(tempIdx) < abs(boundLimits(iB))
                eventBounds(iE,iB) = selectedEvents(iE) - maxBoundRange;
            else
                eventBounds(iE,iB) = candidateBounds(tempIdx(1));
            end
            
            % Positive bound limits, look for the first events after the
            % event time. Take the last one as that is the requested events
            % away.
        elseif boundLimits(iB) > 0
            tempIdx = FastFind(candidateBounds - selectedEvents(iE) > 0,abs(boundLimits(iB)),'first');
            if length(tempIdx) < abs(boundLimits(iB))
                eventBounds(iE,iB) = selectedEvents(iE) + maxBoundRange;;
            else
                eventBounds(iE,iB) = candidateBounds(tempIdx(end));
            end
            
            % Zero bound limits, take the selected behavioral event. Will
            % be after the even for the first limit and before the event
            % for the second limit
        else
            eventBounds(iE,iB) = selectedEvents(iE);
        end
    end
    
    temp = tempTimes < eventBounds(iE,1) | tempTimes > eventBounds(iE,2);
    outOfRange(temp,iE) = 1;  
end

inRange = ~outOfRange;
