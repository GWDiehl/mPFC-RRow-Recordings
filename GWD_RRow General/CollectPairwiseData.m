
function [pairwiseVals, cellPair, cellID] = CollectPairwiseData(fileName,varargin)

dataField = 'Info';
selectedEvent = 'OfferZone';

process_varargin(varargin);


% Collect up all the TE pairwise data
pairwiseData = load(fileName);
Results = pairwiseData.Results;
Params = pairwiseData.Params;

nSess = length(Params.SSN);
pairwiseVals = cell(nSess,1);
cellPair = cell(nSess,1);
cellID = Params.CellID;

% Extract out the normalized TE paired data
cleanPair_Data = cellfun(@(x,y) x-y,Results.(selectedEvent).(dataField),Results.Shuffle.(selectedEvent).Avg,'UniformOutput',0);
for iS = 1:nSess
    nCells = length(cellID{iS});
    CellX = repmat(1:nCells,nCells,1);
    CellY = CellX';
    cellPair{iS} = [CellY(:),CellX(:)];
    pairwiseVals{iS} = cleanPair_Data{iS}(:);
end

% Expand out to full data set as opposed to by session
cellPair = convertSessPairingToFull(cellPair);
pairwiseVals = cell2mat(pairwiseVals);
cellID = cat(1,cellID{1:end});
