function currData = PrepareRRowSingleCellResponses(Results,currDataName,varargin)

% Extracts out the single cell data of interest from a large scale analysis
% structure and normalizes it according to the requested method. Uses the
% ResultID, and OccID fields to extract out the correct naming conventions
% of the data structure.

% GWDiehl 2020


minOcc = .15;
ResultID = 'AvgRate';
OccID = 'AvgOcc';
normMethod = 'None';

% If the data input is 2D+, (cells x D1 x D2 ...) collapse things down into
% a single dimension of data (cells x D_All) to normalize. Otherwise dont
% do any normalization.
collapse2D = 1; 

process_varargin(varargin);

%%


% Reference to the raw data in the case that we are looking at shuffles
genDataName = strrep(currDataName,'Shuff_','');

currData = cat(1,Results.(currDataName).(ResultID){:});
if size(currData,3)>1
    if collapse2D
        fprintf('Your data (%s) has too many dimensions. Collapsing to 1D\n',currDataName)
    else
        fprintf('Skipping %s as it is too many dimensions \n',currDataName)
        return
    end
end

% If we want to exclude by occupancy (passed in field) and the occupancy
% exists (it was calcualted in the first place; if not calculated fair to
% assume no reason to exclude based on occupancy)
if ~isempty(OccID) && isfield(Results.(genDataName),OccID)
    % Expand out Occ according to the number of recorded cells
    fullOcc = expandTCOcc(Results.(genDataName).(OccID),cellfun(@(x) size(x,1),Results.(currDataName).(ResultID)));
    
    % Cat across the full data set
    currOcc = cat(1,fullOcc{:});
    
    % Remove Low Occ bins
    currData(currOcc < minOcc) = nan;
end

% Ensure that data is only 1D
currData = currData(:,:);

switch normMethod
    case 'None'
        % Nothing
    case 'ZScore'
        currData = (currData - nanmean(currData,2))./nanstd(currData,[],2);
    case 'MeanSubtract'
        currData = currData - nanmean(currData,2);
    case 'LogMeanSubtract'
        currData = log(currData);
        currData = currData - nanmean(currData,2);
    case 'ShuffMeanSubtract'
        shuffMean = cat(1,Results.Shuffle.(genDataName).Avg{:});
        currData = currData - shuffMean(:,:);
    case 'LogShuffMeanSubtract'
        shuffMean = cat(1,Results.Shuffle.(genDataName).Avg{:});
        currData = log(currData) - log(shuffMean(:,:));
    case 'ShuffZScore'
        shuffMean = cat(1,Results.Shuffle.(genDataName).Avg{:});
        shuffStd = cat(1,Results.Shuffle.(genDataName).Std{:});
        shuffStd(shuffStd <1e-10) = nan;
        currData = (currData - shuffMean(:,:))./shuffStd(:,:);
end


