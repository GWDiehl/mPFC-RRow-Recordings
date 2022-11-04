

function matchedTimes = MatchRRowEventCount(triggerTimes,subsetParam,LapData,varargin)


matchAcross = 'Site'; % 'Session' or 'Site'

process_varargin(varargin);

%%

if ~iscell(subsetParam)
    subsetParam = {subsetParam};
end
validLaps = IdentifyRRowLaps('All',LapData);
if ~isequal(size(triggerTimes),size(validLaps))
    matchAcross = 'Session';
end

nParams = length(subsetParam);

for iP = 1:nParams    
    tempValid = IdentifyRRowLaps(subsetParam{iP},LapData);
    
    validLaps = validLaps & tempValid;
end

matchedTimes = nan(size(triggerTimes));

switch matchAcross
    case 'Session'
        nEntries = sum(validLaps(:));
        currentEntries = find(~isnan(triggerTimes(:)));
        selectedEntries = datasample(currentEntries,min([nEntries length(currentEntries)]),'replace',false);
        matchedTimes(selectedEntries) = triggerTimes(selectedEntries);
        
    case 'Site'
        for iS = 1:size(validLaps,2)
            nEntries = sum(validLaps(:,iS));
            currentEntries = find(~isnan(triggerTimes(:,iS)));
            selectedEntries = datasample(currentEntries,min([nEntries length(currentEntries)]),'replace',false);
            matchedTimes(selectedEntries,iS) = triggerTimes(selectedEntries,iS);
        end
end





