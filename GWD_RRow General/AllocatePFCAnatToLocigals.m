function groupingIdentity = AllocatePFCAnatToLocigals(FN,cellIDs,anatGroups)

% Takes in a saved file of Anatomy Grouping information and distribute it
% back out into the standard format of logical indexes as would be output
% from the main groupPFCCellsByDepth function

% GWDiehl Oct 2021


anatGrpData = load(FN);

if ~iscell(anatGroups)
    anatGroups = {anatGroups};
end
nCells = length(cellIDs);
nAnat = length(anatGroups);

groupingIdentity = false(nCells,nAnat);

selectedCells = ismember(anatGrpData.CellID,cellIDs);

for iA = 1:nAnat
    if ~isfield(anatGrpData,anatGroups{iA})
        error('Your requested Anatomy group does not exist in the sorted file')
    end
    groupingIdentity(:,iA) = anatGrpData.(anatGroups{iA})(selectedCells);
end