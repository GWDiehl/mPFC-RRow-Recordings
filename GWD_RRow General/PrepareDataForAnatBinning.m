

function [outputResults, Params, outputCounts] = PrepareDataForAnatBinning(pairwiseSim,cellPair,cellID,outputFN,varargin)

% Take a set of pairwise correlation values, and the respective cells and
% prepare the data for binning of the pairwise similairty values according
% to anatomical location. This function largely acts to loop through
% multiple data fields, and subgroup the data ByRat/XRat as needed.
%
% Inputs
% pairwiseSim - Vector of the pairwise similarity values that will be
%       subsequently binned by anatomy. This can be a strucutre with
%       multiple sets of similarity values if you would like to do a set of
%       variables (will have identical anat shuffling across them all).
% cellPair - The Nx2 set of vectors with the cell Idx number for each
%       pairwise value. This point to the respective index in the CellID.
% cellID - A cell array of unique cell names for each cell for which there
%       is a pairing. Index in the cellID array matches the number in the
%       cellPair vector. The total number of entries should match the
%       maximum cellPair entry.
% ouptutFN - What file name would you like output data to be saved under.
%
% Optional
% nBins - How many anat bins do you want (evenly spaced)
% anatBounds - what min/max edges do you want to use
% nShuffles - How many times do you want to shuffle the location data
% shuffMethod - Across what set of data do you want to shuffle location
%   data
% anatGroup - Which set of data do you want to pull anat location data from
% locationIdx - Which index of the location data do you want to look at for
%   binning purposes
% dataSubgrouping - Do you want to subgroup you data in some fashion?
% saveData - logical of save the data to the output Dir?
%
% Outputs
% outputResults - Binned average values
% Params - parameters used for binning the data
% outputCounts - Number of entries in each binned location

% GWD Feb 2021

nBins = 31;
anatBounds = [2.3 5.4];
nShuffles = 30;
shuffMethod = 'Rat';

anatGroup = 'No_VO';
locationIdx = 3; % Which axis of the orig location data do we want to work in
dataSubgrouping = 'All';
symetricData = 1;
groupingMethod = @nanmean;

saveData = 1;

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;

process_varargin(varargin);

%%

% Compute the shuffled depth relationships
[wantedData, stdLocation] = groupPFCCellsByDepth(cellID,anatGroup);
stdLocation(~wantedData,:) = nan; % Remove any entries that are not wanted

shuffLoc = nan([size(stdLocation) nShuffles+1]);
for iX = 1:nShuffles+1
    if iX == 1
        shuffLoc(:,:,iX) = stdLocation;
    else
        shuffLoc(:,:,iX) = shuffleUnitLocation(stdLocation,cellID,'candidates',shuffMethod);
    end
end

% Identify the relevent subgrouping
switch dataSubgrouping
    case {'All' 'None'}
        nSubgroups = 1;
        groupingType = 'ByCell';
        subgrpIdx = true(length(cellID),1);
        
    case 'ByRat'
        unitRat = cellfun(@(x) x(1:4),cellID,'UniformOutput',0);
        subName = unique(unitRat,'stable');
        nSubgroups = length(subName);
        
        groupingType = 'ByCell';
        temp = cellfun(@(x) ismember(unitRat,x),subName,'UniformOutput',0);
        subgrpIdx = cat(2,temp{:});
        
    case {'XRat' 'CrossRat'}
        unitRat = cellfun(@(x) x(1:4),cellID,'UniformOutput',0);
        subName = {'WithinRat' 'CrossRat'};
        nSubgroups = length(subName);
        
        groupingType = 'ByPair';
        sameRat = arrayfun(@(x) isequal(unitRat(cellPair(x,1)),unitRat(cellPair(x,2))),[1:size(cellPair,1)]');
        subgrpIdx = cat(2,sameRat,~sameRat);
        
    case {'XProbe' 'CrossProbe'}
        
        %%% TEST THIS %%%
        
        
        subName = {'WithinProbe' 'CrossProbe'};
        nSubgroups = length(subName);
        
        groupingType = 'ByPair';
        unitRat = cellfun(@(x) x(1:4),cellID,'UniformOutput',0);
        unitProbe = cellfun(@(x) x(end-6:end-3),cellID,'UniformOutput',0);
        sameProbe = arrayfun(@(x) isequal(unitRat(cellPair(x,1)),unitRat(cellPair(x,2))) &&...
            isequal(unitProbe(cellPair(x,1)),unitProbe(cellPair(x,2))),[1:size(cellPair,1)]');
        
        subgrpIdx = cat(2,sameProbe,~sameProbe);
        
    case {'XSess' 'CrossSess'}
        
        %%% TEST THIS %%%
        
        subName = {'WithinSess' 'CrossSess'};
        nSubgroups = length(subName);
        
        groupingType = 'ByPair';
        unitSSN = cellfun(@(x) x(1:15),cellID,'UniformOutput',0);
        sameSSN = arrayfun(@(x) isequal(unitSSN(cellPair(x,1)),unitSSN(cellPair(x,2))),[1:size(cellPair,1)]');
        subgrpIdx = cat(2,sameSSN,~sameSSN);
        
end

% If we have multiple subgroups
if nSubgroups > 1
    FullResults = struct;
    FullCounts = struct;
end

for iG = 1:nSubgroups
    currLoc = shuffLoc;
    
    % Nan out depth of unwanted cells
    if strcmp(groupingType,'ByCell')
        currLoc(~subgrpIdx(:,iG),:) = nan;
    end
    preparedDepth = nan([size(cellPair) nShuffles+1]);
    for iS = 1:nShuffles+1
        for iD = 1:size(cellPair,2)
            preparedDepth(:,iD,iS) = currLoc(cellPair(:,iD),locationIdx,iS);
        end
    end
    % Nan out depth of unwanted pairs
    if strcmp(groupingType,'ByPair')
        preparedDepth(~subgrpIdx(:,iG),:) = nan;
    end
    
    % Remove any entries of preparedDepth in which one of the cells in the
    % pair has been NaNed out (if shuffling was done properly these nan
    % entries should follow across shuffles)
    validEntries = ~any(isnan(preparedDepth(:,:,1)),2);
    preparedDepth = preparedDepth(validEntries,:,:);    
    
    
    if isstruct(pairwiseSim)
        includedField = fieldnames(pairwiseSim);
        nTC = length(includedField);
        Results = struct;
        for iT = 1:nTC
            [tempR, Params, Counts] = binPairwiseAnatSim(pairwiseSim.(includedField{iT})(validEntries),cellPair(validEntries,:),cellID,...
                'nShuffles',nShuffles,'preparedDepth',preparedDepth,'shuffMethod',shuffMethod,...
                'nBins',nBins,'bounds',anatBounds,'symetricData',symetricData,'groupingMethod',groupingMethod);
            
            Results.(includedField{iT}) = tempR;
            fprintf('Done with anatBinning of %s \n',includedField{iT})
        end
    else
        [tempR, Params, Counts] = binPairwiseAnatSim(pairwiseSim(validEntries),cellPair(validEntries,:),cellID,...
            'nShuffles',nShuffles,'preparedDepth',preparedDepth,'shuffMethod',shuffMethod,...
            'nBins',nBins,'bounds',anatBounds,'symetricData',symetricData,'groupingMethod',groupingMethod);
        
        Results.All = tempR;
            fprintf('Done with anatBinning of data \n')
    end
    
    if nSubgroups > 1
        FullResults.(subName{iG}) = Results;
        FullCounts.(subName{iG}) = Counts;
    end
end

if saveData
    if nSubgroups > 1
        save([RootDir,AnalysisDir,outputFN],'FullResults','Params','FullCounts')
    else
        save([RootDir,AnalysisDir,outputFN],'Results','Params','Counts')
    end
end

if nSubgroups > 1    
    outputResults = FullResults;
    outputCounts = FullCounts;
else
    outputResults = Results;
    outputCounts = Counts;
end

fprintf('Done with Pairwise correlations \n')

end


