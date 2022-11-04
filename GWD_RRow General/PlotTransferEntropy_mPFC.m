function [fh,pairwiseTE,cellPair,cellID,outputResults,rawQuantification,outputCounts,binParams] = PlotTransferEntropy_mPFC(TEFileName,varargin)



selectedEventName = 'OfferZone';
cellSubsetSelection = 'All';
logicType = 'AND';

localDist = .2;
minQuantBins = 5;
groupingMethod = @nanmean;
xRange_TE = [-1e-4 6e-4];
xRange_Binned = [-.5e-4 3e-4];

transitionBounds = [];


process_varargin(varargin);

%%

% Collect up all the TE pairwise data
pairwiseData = load(TEFileName);

Results = pairwiseData.Results;
Params = pairwiseData.Params;

nSess = length(Params.SSN);
pairwiseTE = cell(nSess,1);
cellPair = cell(nSess,1);
cellID = Params.CellID;

% Extract out the normalized TE paired data
cleanTE_Data = cellfun(@(x,y) x-y,Results.(selectedEventName).Info,Results.Shuffle.(selectedEventName).Avg,'UniformOutput',0);
for iS = 1:nSess
    nCells = length(cellID{iS});
    CellX = repmat(1:nCells,nCells,1);
    CellY = CellX';
    cellPair{iS} = [CellY(:),CellX(:)];
    pairwiseTE{iS} = cleanTE_Data{iS}(:);
end

% Expand out to full data set as opposed to by session
cellPair = convertSessPairingToFull(cellPair);
pairwiseTE = cell2mat(pairwiseTE);
cellID = cat(1,cellID{1:end});

% Identify the cells to be included in the ensemble

if ~iscell(cellSubsetSelection)
    cellSubsetSelection = {cellSubsetSelection};
end
pairedFields = ismember(cellSubsetSelection,{'SameType' 'DiffType' 'DifferentType'});

if any(pairedFields)
    if ~all(pairedFields)
        validCells = SelectPFCSubset(cellID,cellSubsetSelection(~pairedFields));
    else
        validCells = true(size(cellID));
    end
    switch cellSubsetSelection{pairedFields}
        case 'SameType'
            principal = SelectPFCSubset(cellID,'Principal');
            interneuron = SelectPFCSubset(cellID,'Interneuron');
            
            principalCells(~validCells) = 0;
            interneuron(~validCells) = 0;
            includedPairs = sum(ismember(cellPair,find(principal)),2)==2 | sum(ismember(cellPair,find(interneuron)),2)==2;
            
        case {'DiffType' 'DifferentType'}
            principal = SelectPFCSubset(cellID,'Principal');
            interneuron = SelectPFCSubset(cellID,'Interneuron');
            principal(~validCells) = 0;
            interneuron(~validCells) = 0;
            includedPairs = sum(ismember(cellPair,find(principal)),2)==1 | sum(ismember(cellPair,find(interneuron)),2)==1;
    end
else
    includedCells = SelectPFCSubset(cellID,cellSubsetSelection,'logicType',logicType);
    includedPairs = sum(ismember(cellPair,find(includedCells)),2)==2;
end


% Remove the unwanted TE data
pairwiseTE(~includedPairs) = nan;

% Compute the binned TE Results and then Quantify along the diagonal
[outputResults, binParams, outputCounts] = PrepareDataForAnatBinning(pairwiseTE,cellPair,cellID,[],'symetricData',0,'saveData',0,'nShuffles',0,'groupingMethod',groupingMethod);
pairwiseData = outputResults.All.Raw;
pairwiseData(outputCounts.Raw<minQuantBins) = nan;

nCounts = outputCounts.Raw;
nCounts(nCounts<minQuantBins) = nan;

binCenters = computeBinCenters(binParams.BinEdges{1});
quantBinCenters = linspace(binCenters(1),binCenters(end),length(binCenters)*2-1);
chunkSize = round(localDist/nanmedian(diff(binCenters)));

pairsPerBin = QuantifyMatrixAxis(nCounts,'Diag','chunkSize',chunkSize,'centerMethod',@nansum);

[selectionMean, selectionStd, selectionCount, rawQuantification] = QuantifyMatrixAxis(pairwiseData,'Diag','chunkSize',chunkSize);
selectionMean(selectionCount<minQuantBins) = nan;


% Plot results
SelectDefaultColorMap(SelectCustomColormap('Blue-Yellow-Red'))

fh = figure;
subplot(2,2,1); hold on
histogram(pairwiseTE,linspace(xRange_TE(1),xRange_TE(2),100))
histogram(pairwiseTE,linspace(xRange_TE(1),xRange_TE(2),100),'displaystyle','stairs')
plotVertLine(0,{'--k'})
plotVertLine(nanmean(pairwiseTE),{'r'})
xlim(xRange_TE)
xlabel('Norm Transer Entropy')
ylabel('Number of Pairs')
box off

subplot(2,2,2)
s = imagesc(binCenters,binCenters,pairwiseData);
set(s,'alphadata',~isnan(pairwiseData))
plotIDLine({'color',[.4 .4 .4]})
if isempty(xRange_Binned)
    caxis(prctile(pairwiseData(:),[2.5 97.5]))
else
    caxis(xRange_Binned)
end
colorbar
plotVertLine(transitionBounds,{'--k'})
plotHorizLine(transitionBounds,{'--k'})
xlabel('DV: Input Cell')
ylabel('DV: Receiving Cell')

subplot(2,2,3)
plot(quantBinCenters,selectionMean)
xlim([binParams.BinEdges{1}(1),binParams.BinEdges{1}(end)])
plotVertLine(transitionBounds,{'--k'})
box off
xlabel('Depth (mm)')
ylabel('Normalized TE Strength')


