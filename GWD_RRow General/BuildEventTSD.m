function [outputTSD] = BuildEventTSD(EventTimes,EventData,varargin)

% Creates a standard TSD object off of a set of EventTimes and matching
% EventData. outputTSD will have one entry at each of the EventTimes and a
% corresponding tail entry immedetly before the subsequent EventTime such
% that any utility of TSDs (selecting data etc.) will function normally but
% on a minimally sized TSD. 
%
% Optional Inputs
% stepSize - How far back in time would you like to step from each
%       subsequent EventTime. This should be sufficently small that other
%       events do not fall between the resulting gap as these will behave
%       impropperly (standard problem across all TSDs). Default 1ms.
%
% fullTS - A vector of times or a TS/TSD object with time entries that you
%       would like to have corresponding data for in the outputTSD.
%       Resultant output will match times/data for the input 'fullTS' if
%       you want to fill out your minimal TSD.
% presortTimes - Would you like to ensure that EventTimes are sorted before
%       building the TSD. In genral this should always be the case, though
%       current feature/bugs in TSD constructor code allows for building of out
%       of order TSDs should you have a really good reason doing for such.

% GWD 2019



% Use a stepsize of 1ms. Seems quick and easy and works well for most
% things.
stepSize = 1e-3;

fullTS = [];
presortTimes = 1;

startTime = [];
endTime = [];

process_varargin(varargin);


EventTimes = reshape(EventTimes,[],1);
EventData = reshape(EventData,[],1);

if length(EventTimes) ~=length(EventData)
    error('Times and Data vectors need to be the same length')
end

% Add very first start and end times to include for the TSD that will allow
% for extrapolation outwards from the first/last entry if needed
if ~isempty(startTime)
    EventTimes = [startTime; EventTimes];
    EventData = [EventData(1); EventData];
end
if ~isempty(endTime)
    EventTimes = [EventTimes; endTime];
    EventData = [EventData; EventData(end)];
end

% Presort the event times/data so everything runs sequentially.

% % % In generally this should always be done, however I could imagine
% % % potential instances where it is not wanted. (e.g Previous TSD entry wanted to
% % % be the previous entry in data, not previous timestep.)
if presortTimes
    [EventTimes, sortOrder] = sort(EventTimes);
    EventData = EventData(sortOrder);
end



% Add a buffer entry before so that a preceeding event can be taken to
% bound the first event
EventTimes = [nan; EventTimes];
EventData = [nan; EventData];


% What was the time and event one time step before the current event
PreviousTimeStep = EventTimes(2:end) - stepSize;
PreviousEventData = EventData(1:end-1);

% Group everything and sort
fullTimes = [EventTimes(2:end); PreviousTimeStep];
fullData = [EventData(2:end); PreviousEventData];

% Things have been concatonated to ends so need to resort things
[fullTimes, sortOrder] = sort(fullTimes);
fullData = fullData(sortOrder);


% Make your TSD
outputTSD = tsd(fullTimes,fullData);


% If you want to populate a TSD using different sampling times do that off
% of the partial TSD. Can be a TS, TSD, or simple vector of times

if isempty(fullTS)
    return
elseif isa(fullTS,'ts')
    fullTimes = fullTS.range;
else
    fullTimes = fullTS;
end

if length(fullTimes) < length(outputTSD.range)
    warning('You seem to be further downsampling an already minimal TSD. You will loose data doing this, it is a BAD IDEA.')
end

fullData = outputTSD.data(fullTimes);
outputTSD = tsd(fullTimes,fullData);


