function [requestedData,groupingID,Params] = ExtractPFCAnalysisData(dataFN,dataType,selectedEntry,varargin)

% Loads a set of processed single cell response data for mPFC recordings in
% RRow, identifies the correct field names to extract out the data, calls a
% subsequent function to extract out a selected data field and normalize as
% requested.

% GWDiehl Oct 2021


normMethod = 'ShuffMeanSubtract';

[rootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
[~,minOcc,~,ResultID,OccID] = ResultsFieldNames(dataType);

ensembleData = 0;

process_varargin(varargin);


tempData = load([rootDir,AnalysisDir,dataFN]);
Params = tempData.Params;
Results = tempData.Results;

% If we are wanting shuffled data, pull it up to the top level
if contains(selectedEntry,'Shuff_')
    Results.(selectedEntry).(ResultID) = Results.Shuffle.(strrep(wantedEntries{iW},'Shuff_','')).(ResultID);
end

assert(isfield(Results,selectedEntry),'Your Selected Entry does not seem to exist');


if ensembleData
    groupingID = Params.SSN;
    requestedData = PrepareRRowEnsembleResponse(Results,selectedEntry,...
        'OccID',OccID,'ResultID',ResultID,'normMethod',normMethod,'minOcc',minOcc);    
else  
    groupingID = cat(1,Params.CellID{:});
    requestedData = PrepareRRowSingleCellResponses(Results,selectedEntry,...
        'OccID',OccID,'ResultID',ResultID,'normMethod',normMethod,'minOcc',minOcc); 
    
end

