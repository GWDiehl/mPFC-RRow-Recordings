function mPFC_RevisionAnalysis_GeneralTaskVar_eLife(varargin)


figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\GeneralTask\';
currSubDir = 'FR Relation\';

process_varargin(varargin);

CleanMatlabPlotDefaults

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

%% General Firing rate relation to Maze chunk, Site Rank, and sessionTime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   MazeChunk, Site Rank, Session Time Info   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GWDColor = loadGWD_RRowColormap;
[RootDir, ~, BehavDir, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
dataRootName = 'Revision\eLife\AvgEventRate\';
dataFN = 'GeneralEvents_FullWindow';

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
nShuff = 30;
eventsToUse = {'OZEntry_OZExit','WZEntry_WZExit','LZEntry_LZExit','AZEntry_AZExit'};
spikeBinParams = [5 95 5];
behavBinParams = [5 95 5];

behavVar = {'MazeChunk' 'SessionTime' 'SiteRank'}; %
behavDisc = {[] [0 3600 60] [.5 4.5 4]};
toPlot.Box = true(size(behavVar));
toPlot.Heat = [1 0 0];
toPlot.Line = [1 1 1];

infoRange = [-.05, -.05 -.02 -.05; .8 .45 .17 .35];
anatColorMax = [2.5 2.5 .5 .5; 2.5 2.5 .5 .5; 2.5 2.5 .5 .5];
nBehav = length(behavVar);


fullDataFN = [RootDir,AnalysisDir,dataRootName,dataFN];
temp = load(fullDataFN);
Params = temp.Params;
Results = temp.Results;

nAnat = length(anatGroup);
[subA,subB] = computeSubplotSize(nAnat);
currEvent = eventsToUse;
eventNames = cellfun(@(x) x(1:2),currEvent,'UniformOutput',0);

nEvents = length(currEvent);
prepData = cellfun(@(x) ExtractPFCAnalysisData([dataRootName,dataFN],'General',x,'ensembleData',1),currEvent,'UniformOutput',0);

cellIDs = cat(1,Params.CellID{:});
SSNs = Params.SSN;
nSess = length(SSNs);

% Get any behavioral Task Data that may be needed

[FullLapData, FullSessionData, FullIdPhi] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});
[nSess,nLaps,nSites] = size(FullLapData.Decision);

behavIdxKey = ExtractBehavDataFromFull(FullLapData,[],'LapNumber');

behavData = cell(nBehav,1);
binnedBehav = cell(nBehav,1);
for iB = 1:nBehav
    if ~strcmp(behavVar{iB},'MazeChunk')
        [behavData{iB}, defaultBins] = ExtractBehavDataFromFull(FullLapData,FullSessionData,behavVar{iB});
        
        if isempty(behavDisc{iB})
            binnedBehav{iB} = discretize(behavData{iB},defaultBins);
            behavDisc{iB} = [defaultBins(1) defaultBins(2) length(defaultBins)-1];
        else
            temp = num2cell(behavDisc{iB}); temp{end} = temp{end}+1; % Increment the nBins index to function properly
            binnedBehav{iB} = discretize(behavData{iB},linspace(temp{:}));
        end
    else       
        binnedBehav{iB} = nan(nSess,nLaps,nSites);
        behavDisc{iB} = [.5 nEvents+.5 nEvents];
    end
end

groupID = groupPFCCellsByDepth(cellIDs,anatGroup);
[~,groupIdx] = max(groupID,[],2);
groupIdx(~any(groupID,2)) = nan;
nCells = length(groupIdx);
ratID = cellfun(@(x) x(1:4),cellIDs,'UniformOutput',0);
[allRats,~,ratIdx] = unique(ratID,'stable');
nRats = length(allRats);

% Get the index of the behavioral data that corresponds to each FR event
behavIdx = cell(1,length(currEvent));
behavIdx(:) = {cell(nSess,1);};
for iX = 1:length(currEvent)
    for iS = 1:nSess
        lapIdx = Results.(currEvent{iX}).LapIdx{iS};
        behavIdx{iX}{iS} = arrayfun(@(x) find(behavIdxKey(iS,:)==x),lapIdx);
    end
end
behavBins = cellfun(@(x) x(3),behavDisc);
[behavOutput, avgRate, RawFR, MIData, MI_Shuff] = CollectPerSessFRProfiles(prepData,binnedBehav,behavIdx,behavBins,'nShuff',30);

adjInfo = MIData - MI_Shuff;
adjInfo(isnan(groupIdx),:) = nan;

%% Cross behavior comparison of MI

fh_data = figure;
boxplot(adjInfo,'notch','on','labels',behavVar,'symbol','')
ylim(infoRange(:,end))
ylabel('Norm Info')
box off

figFN = [figOutputDir,currSubDir,sprintf('MI_AcrossMetrics')];
saveas(fh_data,figFN,'svg'); close(fh_data)

allP = nan(size(adjInfo,2),1);
allStats = cell(size(adjInfo,2),1);
for iX = 1:size(adjInfo,2)
    [allP(iX),h,allStats{iX}] = signrank(adjInfo(:,iX));
end

[p,ANOVA,ANOVA_stats] = anova1(adjInfo,[],'off');
multiCompare = multcompare(ANOVA_stats,'display','off');

%% Boxplot and error bar of MI across mPFC subregions

allP = nan(nBehav,nAnat);
allStats = cell(nBehav,nAnat);
for iB = 1:nBehav
    if toPlot.Box(iB);
        currInfo = adjInfo(:,iB);
        fh_data = figure; hold on
        boxplot(currInfo,groupIdx,'notch','on','labels',anatGroup,'symbol','.','jitter',.3)
        
        [p,ANOVA,ANOVA_stats] = anova1(currInfo,groupIdx,'off');
        multiCompare = multcompare(ANOVA_stats,'display','off');
        
        meanData = nan(nAnat,nRats);
        for iA = 1:nAnat
            [allP(iB,iA),~,allStats{iB,iA}] = signrank(currInfo(groupID(:,iA)));
            for iR = 1:nRats
                selectIdx = groupID(:,iA) & ratIdx == iR;
                meanData(iA,iR) = nanmean(currInfo(selectIdx));
            end
        end
        errorbar(1:nAnat,nanmean(meanData,2),nanstderr(meanData'),'k','linewidth',1);
        ylim(infoRange(:,iB))
        box off
        ylabel(sprintf('Norm MI to %s',behavVar{iB}))
        
        figFN = [figOutputDir,currSubDir,sprintf('MI_%s',behavVar{iB})];
        saveas(fh_data,figFN,'svg'); close(fh_data)
    end
end

%% Color scale/Line Plot of the avg rate in each behavior variable bin

for iB = 1:nBehav    
    currFR = avgRate{iB};
    meanData = nan(nAnat,behavBins(iB));
    errorData = nan(nAnat,behavBins(iB));
    for iA = 1:nAnat
        for iE = 1:behavBins(iB)
            meanData(iA,iE) = nanmean(currFR(groupID(:,iA),iE));
            errorData(iA,iE) = nanstderr(currFR(groupID(:,iA),iE));
        end
    end
    if toPlot.Heat(iB)
        fh_data = figure;
        for iA = 1:nAnat
            subplot(subA,subB,iA)
            imagesc(meanData(iA,:))
            if strcmp(behavVar{iB},'MazeChunk')
                xticks([1:nEvents]);xticklabels(eventNames)
            end
            if ~isnan(anatColorMax(iB,iA));
                caxis([-anatColorMax(iB,iA) anatColorMax(iB,iA)])
            end
            xlabel(behavVar{iB});
            title(anatGroup{iA})
            colorbar
        end
        figFN = [figOutputDir,currSubDir,sprintf('AvgFR_%s_Heatmap',behavVar{iB})];
        saveas(fh_data,figFN,'svg'); close(fh_data)
    end
    if toPlot.Line(iB)
        fh_data = figure; hold on
        temp = num2cell(behavDisc{iB}); temp{end} = temp{end}+1;
        for iA = 1:nAnat
            errorbar(computeBinCenters(linspace(temp{:})),meanData(iA,:),errorData(iA,:),'color',GWDColor.(anatGroup{iA}));
            plot(computeBinCenters(linspace(temp{:})),meanData(iA,:),'o','color',GWDColor.(anatGroup{iA}))
        end
        xlim([temp{1},temp{2}]);
        xlabel(behavVar{iB})
        ylabel('Norm FR')
        if strcmp(behavVar{iB},'MazeChunk')
            xticks([1:nEvents]);xticklabels(eventNames)
        end
        figFN = [figOutputDir,currSubDir,sprintf('AvgFR_%s_LinePlot',behavVar{iB})];
        saveas(fh_data,figFN,'svg'); close(fh_data)
    end
end

