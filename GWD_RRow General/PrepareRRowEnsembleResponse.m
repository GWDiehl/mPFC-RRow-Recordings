function currData = PrepareRRowEnsembleResponse(Results,currDataName,varargin)

% Extracts out the single cell data of interest from a large scale analysis
% structure and normalizes it according to the requested method. Uses the
% ResultID, and OccID fields to extract out the correct naming conventions
% of the data structure.

% GWDiehl Oct 2021

minOcc = .15;
ResultID = 'AvgCorr';
OccID = 'AvgOcc';
normMethod = 'None';

process_varargin(varargin);

%%


% Reference to the raw data in the case that we are looking at shuffles
genDataName = strrep(currDataName,'Shuff_','');

currData = Results.(currDataName).(ResultID);

for iX = 1:numel(currData)
    % Remove Low Occ bins
    % If we want to exclude by occupancy (passed in field) and the occupancy
    % exists (it was calcualted in the first place; if not calculated fair to
    % assume no reason to exclude based on occupancy)
    if ~isempty(OccID) && isfield(Results.(genDataName),OccID)
        fullOcc = Results.(genDataName).(OccID);
        currData{iX}(fullOcc{iX}<minOcc,:) = nan;        
        currData{iX}(:,fullOcc{iX}<minOcc) = nan;
    end
    
    % Normalize each entry as requested
    switch normMethod
        case 'None'
            % Nothing
        case 'ZScore'
            currData{iX} = (currData{iX} - nanmean(currData{iX}(:)))./nanstd(currData{iX}(:));
        case 'MeanSubtract'
            currData{iX} = currData{iX} - nanmean(currData{iX},2);
        case 'LogMeanSubtract'
            currData{iX} = log(currData{iX});
            currData{iX} = currData{iX} - nanmean(currData{iX}(:));
        case 'ShuffMeanSubtract'
            shuffMean = Results.Shuffle.(genDataName).Avg{iX};
            currData{iX} = currData{iX} - shuffMean;
        case 'LogShuffMeanSubtract'
            shuffMean = Results.Shuffle.(genDataName).Avg{iX};
            currData{iX} = log(currData{iX}) - log(shuffMean);
        case 'ShuffZScore'
            shuffMean = Results.Shuffle.(genDataName).Avg{iX};
            shuffStd = Results.Shuffle.(genDataName).Std{iX};
            shuffStd(shuffStd <1e-10) = nan;
            currData{iX} = (currData{iX} - shuffMean)./shuffStd;
    end
end

