function selectedCells = SelectPFCSubset(cellIDs,selectedGroup,varargin)


% Receives a set of cellIDs in the standard name format and identifies a
% specific subset of those cells based on a variety of criteria entered in
% 'selectedGroup'. SelectedGroup should be a cell array of entries and can
% accomidate multiple different entry fields. 
%
% Defaults to useing  AND logic between them, but can also use OR logic. 
% **Double check that exclusions option works properly with the OR logic**

% GWDiehl Oct 2021

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] ...
    = ID_GWD_Dirs;

cellTypeFN = [RootDir,AnalysisDir,'\CellTypeClassification'];
anatRegionFN = [];

logicType = 'AND';  % AND  OR

process_varargin(varargin);

%%

if ~iscell(selectedGroup)
    selectedGroup = {selectedGroup};
end

cellTypeInfo = load(cellTypeFN);

nCells = length(cellIDs);
switch logicType
    case 'AND'
        selectedCells = true(nCells,1);
    case 'OR'
        selectedCells = false(nCells,1);
end

nGroupings = length(selectedGroup);

ratID = cellfun(@(x) x(1:4),cellIDs,'UniformOutput',0);
sessID = cellfun(@(x) x(1:15),cellIDs,'UniformOutput',0);

for iG = 1:nGroupings
    
    currCriteria = selectedGroup{iG};
    excludeData = 0;
    
    if strcmp(currCriteria(1),'~')
        excludeData = 1;
        currCriteria = currCriteria(2:end);
    end    
    
    switch currCriteria
        case 'All'
            currSelection = true(size(cellIDs));
        case {'No_VO' 'Acc' 'Acc_TE' 'PL' 'dPL' 'dPL_V2' 'dPL_TE' 'vPL' 'vPL_V2' 'vPL_TE' 'IL' 'IL_TE'}
            if isempty(anatRegionFN)
                currSelection = groupPFCCellsByDepth(cellIDs,currCriteria);
            else % If there exists a file that keys to which cell was in what area, just use that
                currSelection = AllocatePFCAnatToLocigals(anatRegionFN,cellIDs,currCriteria);
            end                
        case {'Principal'}
            cellSubset = ismember(cellTypeInfo.CellID,cellIDs);
            currSelection = cellTypeInfo.Principal(cellSubset);
        case {'Interneuron'}
            cellSubset = ismember(cellTypeInfo.CellID,cellIDs);
            currSelection = cellTypeInfo.Interneuron(cellSubset);
            
        otherwise
            % Current Criteria is a Rat
            if length(currCriteria)==4 && strcmp(currCriteria(1),'R')
                currSelection = ismember(ratID,currCriteria);
                % Session ID (SSN)
            elseif length(currCriteria) == 15
                currSelection = ismember(sessID,currCriteria);
                % Full Cell ID
            elseif length(currCriteria) == 23
                currSelection = ismember(cellIDs,currCriteria);
            else
                error('Case Not Specified')
            end
    end
    
    if excludeData
        currSelection = ~currSelection;
    end
    
    switch logicType
        case 'AND'            
            selectedCells = selectedCells & currSelection;
        case 'OR'
            selectedCells = selectedCells | currSelection;
    end
end
