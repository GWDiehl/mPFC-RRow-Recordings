
function [fhScatter,fhBox,fullData] = BasicSpikingVsAnat_V2(allSpks,cellID,spikingQuants,varargin)


% Compute and plot basic spiking characteristics of recorded mPFC cells


% GWD Nov 2021

sessDuration = 60*60;
minSpks = 50;
depthBin = 3;

nBins = 31;
anatBounds = [2.3 5.4];
depthRange = [2.3 5.4];
minCells = 5;

cellGroupings = {'No_VO'};
groupNames = {'No_VO'};
plotScatter = [1];
plotErrorbar = [1];
GWDColorMap = loadGWD_RRowColormap;

process_varargin(varargin);

%%

nGrps = length(cellGroupings);
anatEdges = linspace(anatBounds(1),anatBounds(2),nBins+1);

nCells = length(allSpks);
nGrps = length(cellGroupings);
nQuants = length(spikingQuants);

if any(ismember(spikingQuants,{'PeakValleyRatio' 'SpikeWidth'}))
    cellWF = ExtractCellWF_GWD(cellID);
else
    cellWF = [];
end

% Find the depth info and remove anything outside of the bounds of interest
stdLocation = FindUnitLocation_mPFC(cellID);
stdDepth = stdLocation(:,depthBin);
stdDepth(stdDepth < depthRange(1) | stdDepth > depthRange(2)) = nan;

fhScatter = cell(nQuants,1);
fhBox = cell(nQuants,1);
fullData = nan(nCells,nQuants,nGrps);

nSpks = cellfun(@(x) length(x.range),allSpks);

selectedCells = false(nCells,nGrps);
for iG = 1:nGrps
    selectedCells(:,iG) = SelectPFCSubset(cellID,cellGroupings{iG});
end

for iQ = 1:nQuants
    
    [basicQuant,defaultRange] = QuantifyBasicCellSpiking(allSpks,spikingQuants{iQ},cellWF);
    basicQuant(nSpks<minSpks) = nan;
    
    fhScatter{iQ} = figure;hold on
    tempData = [];
    tempGrp = [];
    
    for iG = 1:nGrps
        
        currSelection = selectedCells(:,iG);
        
        currDepth = stdDepth(currSelection);
        currData = basicQuant(currSelection);
        valid = ~isnan(currDepth) & ~isnan(currData);
        
        if plotScatter(iG)
            if ~iscell(cellGroupings{iG}) && isfield(GWDColorMap,cellGroupings{iG})
                scatter(currDepth(valid),currData(valid),'MarkerEdgeColor',GWDColorMap.(cellGroupings{iG}))
            else
                scatter(currDepth(valid),currData(valid))
            end
        end
        
        [dataMean, count] = BinPairWiseData(currData(valid),currDepth(valid),'Bins',{anatEdges},'groupingMethod',@nanmean);
        [dataSEM] = BinPairWiseData(currData(valid),currDepth(valid),'Bins',{anatEdges},'groupingMethod',@nanstderr);
        dataMean(count<minCells) = nan;
        
        if plotErrorbar(iG)
            errorbar(computeBinCenters(anatEdges),dataMean,dataSEM,'k','LineWidth',1)
        end
        
        fullData(currSelection,iQ,iG) = currData;
        tempData = cat(1,tempData,currData(valid));
        tempGrp = cat(1,tempGrp,repmat(iG,sum(valid),1));
    end
    
    xlabel('Std Depth')
    ylabel(spikingQuants{iQ})
    ylim(defaultRange)
    
    
    fhBox{iQ} = figure;
    boxplot(tempData,tempGrp,'Labels',groupNames,'notch','on')
    ylabel(spikingQuants{iQ})
    box off
    ylim(defaultRange)
end


