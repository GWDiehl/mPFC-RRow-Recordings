function mPFC_RevisionAnalysis_TransferEntropy_eLife(varargin)


figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\';

process_varargin(varargin);


CleanMatlabPlotDefaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Correlation to specific behavioral vars   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'TransferEntropy\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

%%

[rootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;

inputBase = 'Revision\eLife\TransferEntropy\MazeChunks\';
TEFileName = 'PairwiseTE_FullEvent';
mazeChunk = {'All' 'OfferZone' 'WaitZone' 'LingerZone' 'AdvanceZone'};

localDist = .2;
transitionBounds = [2.9 3.85 4.45];
cellSubsetSelection = 'No_VO';
minQuantBins = 5;

tempResultsOutput = cell(length(mazeChunk),1);

fh_full = figure; hold on
fullOutput = nan(31,31,length(mazeChunk));
for iX = 1:length(mazeChunk)
    fullFN = [rootDir,AnalysisDir,inputBase,TEFileName,'_',mazeChunk{iX}];
    [fh,pairwiseTE,cellPair,cellID,outputResults,rawQuantification,outputCounts,binParams] = ...
        PlotTransferEntropy_mPFC(fullFN,'cellSubsetSelection',cellSubsetSelection,'transitionBounds',transitionBounds,'selectedEventName',mazeChunk{iX});
    
    binnedData = outputResults.All.Raw;
    binnedData(outputCounts.Raw<5) = nan;
     fullOutput(:,:,iX) = binnedData;
    
    [r, p] = nancorr(binnedData,binnedData');
    
    [p, ~, stats] = signrank(pairwiseTE);
    
    binCenters = computeBinCenters(binParams.BinEdges{1});
    quantBinCenters = linspace(binCenters(1),binCenters(end),length(binCenters)*2-1);
    
    selectionMean = nanmean(rawQuantification,2);
    selectionCount = sum(~isnan(rawQuantification),2);
    selectionMean(selectionCount<minQuantBins) = nan;
    [pks, pkIdx,w,p] = findpeaks(-selectionMean,quantBinCenters,'MinPeakProminence',5e-5);
    
    % Run a Hartigan Dip Test for unimodal significance on the quantificaition
    validBins = ~isnan(selectionMean);
    sampleData = randsample(quantBinCenters(validBins),1000,true,selectionMean(validBins));
    [Dip, p_value]=HartigansDipSignifTest(sampleData,500);
    
    figFN = [figOutputDir,currSubDir,'TE_ByDepth_',mazeChunk{iX}];
    fh.Position = [100 100 1150 620];
    saveas(fh,figFN,'svg'); close(fh)
    
    figure(fh_full)
    plot(quantBinCenters,selectionMean)
    tempResultsOutput{iX} = selectionMean;
end

box off
xlabel('Depth (mm)')
ylabel('Normalized TE Strength')
xlim([binParams.BinEdges{1}(1),binParams.BinEdges{1}(end)])
plotVertLine(transitionBounds,{'--k'})
legend(mazeChunk,'location','best')

figFN = [figOutputDir,currSubDir,'TE_ByDepth_AllChunks'];
fh_full.Position = [100 100 1150 620];
saveas(fh_full,figFN,'svg'); close(fh_full)

compOutput = nan(length(mazeChunk));
for iX = 1:length(mazeChunk)
    for iY = 1:length(mazeChunk)
        compOutput(iX,iY) = nancorr(fullOutput(:,:,iX),fullOutput(:,:,iY));
    end
end
minCorr = min(compOutput(:));


%% Compute TE distributions within a subregion and between different combinations

[rootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;

inputBase = 'Revision\eLife\TransferEntropy\MazeChunks\';
TEFileName = 'PairwiseTE_FullEvent';
mazeChunk = 'All';

fullFN = [rootDir,AnalysisDir,inputBase,TEFileName,'_',mazeChunk];
anatGroups = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};

[pairwiseVals, cellPair, cellID] = CollectPairwiseData(fullFN,'selectedEvent',mazeChunk);

[groupedTEVals,fh_byRegion,fh_byDistance] ...
    = PlotTEAcrossSubregions(pairwiseVals,cellPair,cellID,anatGroups);

% figFN = [figOutputDir,currSubDir,'TE_AcrossSubregions'];
% fh_byRegion.Position = [100 100 1150 620];
% saveas(fh_byRegion,figFN,'svg'); 
close(fh_byRegion)

figFN = [figOutputDir,currSubDir,'TE_WithinVsAcross'];
fh_byDistance.Position = [100 100 1150 620];
saveas(fh_byDistance,figFN,'svg'); close(fh_byDistance)


%% DV Location info of the atlas grouped Acc, PL, IL cells

anatGroup = {'Acc' 'PL' 'IL'};
anatBounds = [2.3 5.4];
transitionBounds = [2.9 3.85 4.45];
nBins = 31;
anatEdges = linspace(anatBounds(1),anatBounds(2),nBins+1);
anatCenters = computeBinCenters(anatEdges);

grpColors = loadGWD_RRowColormap;

[~,~,~,~, Full_FNs] = GWD_LoadAllRRowSessions('dataToLoad',{'SpikeNames'});
S_FNs = cat(1,Full_FNs{:});

[dataGroup,stdLocation] = groupPFCCellsByDepth(S_FNs,anatGroup);

fh = figure;hold on
for iG = 1:length(anatGroup)
    cellCount = histcn(stdLocation(dataGroup(:,iG),3),anatEdges);
    plot(anatCenters,cellCount,'color',grpColors.(anatGroup{iG}))
end
xlim(anatBounds);
xlabel('DV Location')
ylabel('Cell Count')
title('Grouped by Atlas Designations')

figure(fh); plotVertLine(transitionBounds,{'--k'})
legend(anatGroup)
figFN = [figOutputDir,currSubDir,'Atlas_DVLoc'];
saveas(fh,figFN,'svg'); close(fh)


