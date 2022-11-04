function shuffledLoc = shuffleUnitLocation(origLocation,cellID,varargin)

% Take a set of 3D single unit locations and shuffle them among themselves
% according to one of a set of different methods: All Data, Within the Same
% Probe, Within the same Rat, etc.
%
% Inputs
% origLocation - nCells x 3 mtx of the location of each cell.
% celID - Cell array of the corresponding cell identifier composed of
%       SSN-TTXX_XX or SSN-SiXX_XX. Note: All data will behave
%       appropriately so long as it has a 15 char SSN followed by the above
%       Seneor/Unit designation.
%
% Outputs
% shuffledLoc - A rearanged set of 3D location data.
%
% Optional
% removeUnshuffled - NaN out any data points that lacked an appropriate
%       shuffle candidate from the data set. Alternatively these entries
%       can be returned unchanged.
% candidates - Str of the appropriate candidate pool you would like to pull
%       from. ('All' 'Probe' 'Rat' 'Session' 'Session+Probe')

% GWD May 2020


removeUnshuffled = 1;
candidates = 'All';

process_varargin(varargin);


if iscell(origLocation(1))
    breakCell = 1;
    fullCellID = cat(1,cellID{:});
    fullLocation = cat(1,origLocation{:});
else
    breakCell = 0;
    fullLocation = origLocation;
    fullCellID = cellID;
end

% Get the info you need
switch candidates
    case 'Probe' % Any recording from that probe
        unitProbe = cellfun(@(x) x(end-6:end-3),fullCellID,'UniformOutput',0);
        unitRat = cellfun(@(x) x(1:4),fullCellID,'UniformOutput',0);
        [~,~,unitProbe] = unique(unitProbe,'stable');
        [~,~,unitRat] = unique(unitRat,'stable');
        
    case 'Rat' % Any recording from that rat
        unitRat = cellfun(@(x) x(1:4),fullCellID,'UniformOutput',0);
        [~,~,unitRat] = unique(unitRat,'stable');
        
    case 'Session' % Any recording from that session
        unitSSN = cellfun(@(x) x(1:15),fullCellID,'UniformOutput',0);
        [~,~,unitSSN] = unique(unitSSN,'stable');
        
    case {'Session+Probe' 'Probe+Session'} % Any recording from that probe in that session
        unitProbe = cellfun(@(x) x(end-6:end-3),fullCellID,'UniformOutput',0);
        unitSSN = cellfun(@(x) x(1:15),fullCellID,'UniformOutput',0);
        [~,~,unitProbe] = unique(unitProbe,'stable');
        [~,~,unitSSN] = unique(unitSSN,'stable');
end

nCells = length(fullCellID);
newIdx = nan(nCells,1);

voidCells = isnan(fullLocation(:,1));
fullLocation(end+1,:) = nan;
unitNum = 1:nCells;

for iC = 1:nCells
    
    if voidCells(iC) % If the particular cell is voided disregard it
        candidate = 0;
    else        
        switch candidates
            case 'All' % Anything recoring the input dataset
                candidate = true(size(newIdx));
                
            case 'Probe' % Any recording from that probe
                candidate = (unitProbe == unitProbe(iC)) & (unitRat == unitRat(iC));
                
            case 'Rat' % Any recording from that rat
                candidate = (unitRat == unitRat(iC));
                
            case 'Session' % Any recording from that session
                candidate = (unitSSN == unitSSN(iC));
                
            case {'Session+Probe' 'Probe+Session'} % Any recording from that probe in that session
                candidate = (unitProbe == unitProbe(iC)) & (unitSSN == unitSSN(iC));
            otherwise
                error('Please use a valid shuffle method');
        end
        
        candidate(iC) = false; % Cant be the same cell
        candidate(voidCells) = false; % Dont select a voided cell
    end
    
    if sum(candidate) == 0
        if removeUnshuffled % No other candidates, just disregard the cell (use Nan)
            newIdx(iC) = nCells+1;
        else % Use its unshuffled location
            newIdx(iC) = iC;
        end
    else
        newIdx(iC) = randsample(unitNum(candidate),1);
    end
    
end

tempLocation = fullLocation(newIdx,:);

if breakCell
    nGroups = length(cellID);
    shuffledLoc = cell(nGroups,1);
    for iG = 1:nGroups
        count = length(cellID{iG});
        shuffledLoc{iG} = tempLocation(1:count,:);
        tempLocation(1:count,:) = [];
    end
else
    shuffledLoc = tempLocation;
end

