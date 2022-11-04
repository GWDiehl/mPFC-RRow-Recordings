function [Results, Params, Counts] = binPairwiseAnatSim(similarityData,pairing,cellID,varargin)

% Take a set of single metric comparisions between a pair (or more/less) of
% cells and avg/bin those data according to the anatomical depth of the
% pair of cells involved. Originally built for Pairwise correlation between
% FR (QMtx) or TuningCurves.
%
% INPUTS
% similairtyData - A N x 1 vector of all of the pairwise comparisons
% pairing - An N x 2(+ or -) Mtx of the corresponding unit identifiers for
%       those cells involved in the similairty metric.
% cellID - A list of the candidate units where the number in 'pairing'
%       corresponds to the indx of the respective cell in 'cellID'
% NOTE: All three of these can be also input as cell arrays in the case
%       where calculations are done for each session independently (i.e.
%       reseting set of pairing idx each session)
%
% OUTPUT
% Results - A structure containing the desired binned histograms
%       (and shuffles) according to if you want varients of the similairty
%       data (Abs, Pos, etc.)
% Counts - Corresponding structure to Results with the number of entries in
%       each hist bin. Note that Results from raw data are complete and not
%       excluded based on count thresholds (do this posthoc). Shuffles can
%       based on logical.
% Params - General parameters of the analyses.
%
% Optional
% cellsToExclude - A logical array corresponding to 'cellID' of any entries
%       that you would like excluded from the calculation. These will not
%       factor into hist calc or to shuffling of anat location (Nan out).
% dataToExclude - A logical array corresponding to 'similarityData' &
%       'pairing' of any entries that you would like excluded. These will
%       be removed from the data and not included in hist/counts, but will
%       not have any impact of shuffling of anat location.
% nBins - Number of bins bounds - Boundary conditions. Both nBins and
%       bounds should match number of cells in pairing (size,2)
% symetricData - If the axis are equivelent collapose everything to one
%       side of the diag
% polaritySplit - Cell array of which sets of polarity splits of the
% similarity data you would
%       like to compute (Raw, Abs, Pos, Neg)
% nShuffles - How many Shuffles shuffMethod - How do you want to be
%       shuffling the anat data.
% keepAllShuffles - Keep the full set of shuffled data or just Mean/Std
% removeShuffUndersample - Remove any bins without enough samples before
%       Mean/Std of shuffles.
% preparedDepth - A user input set of precalculated depth metrics. In the
%       event of running a series of comparisons with identical shuffling
%       procedures this approach precludes the need to redo the shuffles
%       each time which can be long. NOTE: Input 'preparedDepth' should
%       always match 'similarityData' & 'pairing'.


% GWD Oct 2020


%%

cellsToExclude = [];
dataToExclude = [];

nBins = [31];
bounds = [2.3 5.4];

% Are data symetric about the diagonal equivelent? In that case we will
% collapse them all to one side
symetricData = 1;

polaritySplit = {'Raw' 'Abs' 'Positive' 'Negative'};

nShuffles = 0;
shuffMethod = 'Rat'; %
removeShuffUndersample = 1;
minCounts = 5;
keepAllShuffles = 0;
fullPairedDepth = [];
fullShuffDepth = [];

groupingMethod = @nanmean;

process_varargin(varargin);


%% Binning for making histograms
nBins = nBins(:);
if numel(nBins) == 1
    nBins = repmat(nBins,2,1);
end
if size(bounds,1) == 1
    bounds = repmat(bounds,2,1);
end
Bins = cell(1,length(nBins));
for iD = 1:length(nBins)
    Bins{iD} = linspace(bounds(iD,1), bounds(iD,2), nBins(iD)+1);
end

%% Initialize Results

% What are we doing as far as the depth information? Really only relevant
% regarding shuffles becuase otherwise it is what it is.
if ~isempty(fullPairedDepth)
    % Provided depth info for all cell pairings
    depthSetting = 'PairsProvided';
    assert(size(fullPariedDepth,3)==nShuffles+1,'You did not input enough depth data')
    
elseif ~isempty(fullShuffDepth)
    % Provided depth info for each cell but need to be joined to the
    % pairings
    depthSetting = 'CellsProvided';    
    assert(size(fullShuffDepth,3)==nShuffles+1,'You did not input enough depth data')
    
else
    % Compute the depth info internally
    depthSetting = 'Computed';
end


Results = struct;
Counts = struct;

Params = struct;
Params.nBins = nBins;
Params.Bounds = bounds;
Params.BinEdges = Bins;
Params.nShuffles = nShuffles;

Params.DepthInfo = depthSetting;
Params.CellExclusion = cellsToExclude;
Params.DataExclusion = dataToExclude;
if nShuffles
    Params.RemoveShuffUndersample = removeShuffUndersample;
    Params.minCounts = minCounts;
    Params.ShuffleMethod = shuffMethod;
    Params.KeepAllShuffles = keepAllShuffles;
end

%% Prep data as needed

% If data were originally calculated per session and come in as a cell
% array, adapt it to be a full data set formating.
if iscell(similarityData(1))
    pairing = convertSessPairingToFull(pairing);
    similarityData = cell2mat(similarityData);
    cellID = cat(1,cellID{:});
    if ~isempty(cellsToExclude)
        cellsToExclude = cat(1,cellsToExclude{:});
    end
    if ~isempty(dataToExclude)
        dataToExclude = cat(1,dataToExclude{:});
    end
end

% Remove (full remove) any unwanted paired data
if ~isempty(dataToExclude)
    similarityData = similaritydata(~dataToExclude);
    pairing = pairing(~dataToExclude,:);
end

%% Gather the respective depth info

% Any shuffling is done inside this fuction
switch depthSetting
    case 'PairsProvided'
        if ~isempty(dataToExclude)
            fullPairedDepth = fullPairedDepth(~dataToExclude,:,:);
        end
        
    case 'CellsProvided'
        if ~isempty(cellsToExclude)
            fullShuffDepth(cellsToExclude,:) = nan;
        end
        
        fullPairedDepth = nan(size(pairing));
        for iD = 1:size(pairing,2)
            fullPairedDepth(:,iD) = fullShuffDepth(pairing(:,iD,1),3);
        end
        
        % Unused by required for parfor to function
        stdLocation = [];
        
    case 'Computed'
        stdLocation = FindUnitLocation_mPFC(cellID);
        % Remove (nanout) any unwanted cells
        if ~isempty(cellsToExclude)
            stdLocation(cellsToExclude,:) = nan;
        end
        
        fullPairedDepth = nan(size(pairing));
        for iD = 1:size(pairing,2)
            fullPairedDepth(:,iD) = stdLocation(pairing(:,iD),3);
        end
end

if symetricData % Ordering of dimensions does not matter
    fullPairedDepth = sort(fullPairedDepth,2);
end


%% Loop through and make the histograms

for iP = 1:length(polaritySplit)
    dataIn = similarityData;
    switch polaritySplit{iP}
        case 'Raw'
            % All Good
            validIdx = true(size(dataIn,1),1);
        case 'Abs'
            dataIn = abs(dataIn);
            validIdx = true(size(dataIn,1),1);
        case 'Positive'
            validIdx = dataIn>0;
        case 'Negative'
            validIdx = dataIn<0;
        otherwise
            error('Polarity case not defined')
    end
    
    [Results.(polaritySplit{iP}), Counts.(polaritySplit{iP})] = BinPairWiseData(dataIn(validIdx),fullPairedDepth(validIdx,:,1),'Bins',Bins,'groupingMethod',groupingMethod);
    
    
    if nShuffles
        tempShuff = nan([nShuffles nBins']);
        
        switch depthSetting
            case 'PairsProvided'
                for iS = 1:nShuffles
                    shuffDepth = fullPairedDepth(:,:,iS+1);
                    [binnedData, pairsCount] = BinPairWiseData(dataIn(validIdx),fullPairedDepth(validIdx,:,iS+1),'Bins',Bins,'groupingMethod',groupingMethod);
                    if removeShuffUndersample
                        binnedData(pairsCount < minCounts) = nan;
                    end
                    tempShuff(iS,:) = binnedData(:);
                end
                
            case {'CellsProvided' 'Computed'}  
                parfor iS = 1:nShuffles                    
                    % Get the particular depth info (passed in or shuffled
                    % here)                    
                    tempLoc = [];
                    switch depthSetting
                        case 'CellsProvided'
                            tempLoc = fullShuffDepth(:,:,iS+1);
                        case 'Computed'
                            tempLoc = shuffleUnitLocation(stdLocation,cellID,'candidates',shuffMethod);
                    end
                    
                    % Go do the stuff
                    shuffDepth = nan(size(pairing));                    
                    for iD = 1:size(pairing,2)
                        shuffDepth(:,iD) = tempLoc(pairing(:,iD),3);
                    end
                    if symetricData % Ordering of dimensions does not matter
                        shuffDepth = sort(shuffDepth,2);
                    end
                    [binnedData, pairsCount] = BinPairWiseData(dataIn(validIdx),shuffDepth(validIdx,:),'Bins',Bins,'groupingMethod',groupingMethod);
                    if removeShuffUndersample
                        binnedData(pairsCount < minCounts) = nan;
                    end
                    tempShuff(iS,:) = binnedData(:);
                end
        end
        
        if keepAllShuffles
            Results.Shuffle.(polaritySplit{iP}) = tempShuff;
        else
            Results.Shuffle.(polaritySplit{iP}).Avg = reshape(nanmean(tempShuff,1),nBins');
            Results.Shuffle.(polaritySplit{iP}).Std = reshape(nanstd(tempShuff,1),nBins');
        end
    end
end





