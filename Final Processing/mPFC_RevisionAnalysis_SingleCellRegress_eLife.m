function mPFC_RevisionAnalysis_SingleCellRegress_eLife(varargin)

figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\';

process_varargin(varargin);


CleanMatlabPlotDefaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Correlation to specific behavioral vars   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'SingleCellRegression\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

%%

% Load in the stepwise regression data and get cellIDs
[RootDir, ~, BehavDir, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
GWDColor = loadGWD_RRowColormap;
setJetAsDefaultColorMap

dataDirRoot = 'Revision\eLife\SingleCellRegression\';
dataFN = 'BasicRelationToBehavior';

% Identify which subregion each cell comes from
anatGroup = {'No_VO' 'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
fieldsToExclude = {'Random' 'Shuffle'};

temp = load([RootDir,AnalysisDir,dataDirRoot,dataFN]);

Results = temp.Results;
Params = temp.Params;

varOfInterest = Params.BehavVar;
for iX = 1:length(fieldsToExclude)
    varOfInterest(ismember(varOfInterest,fieldsToExclude{iX})) = [];
end
nVar = length(varOfInterest);

allColors = jet(nVar);
nSig = nan(nVar);

dataCorr = ExtractRegressionOutputs(Results,'Individual','AbsCorrelation');
data_InModel = ExtractRegressionOutputs(Results,'Individual','InModel');
shuffResults = SimplifiyStructByEntries(Results.Shuffle,'Example');
shuffCorr = ExtractRegressionOutputs(shuffResults,'Individual','AbsCorrelation');

% Do each anatGroup seperatly
for iA = 1:length(anatGroup)
    groupID = groupPFCCellsByDepth(cat(1,Params.CellID{:}),anatGroup(iA));
    
    currSubDir = ['SingleCellRegression\',anatGroup{iA},'\'];
    
    if ~exist([figOutputDir,currSubDir],'dir')
        mkdir([figOutputDir,currSubDir]);
    end
    
    fh_data = figure; hold on    
    statsResults = cell(nVar,1);
    pResults = nan(nVar,1);
    for iV = 1:nVar
        behavIdx = ismember(Params.BehavVar,varOfInterest{iV});
        currData = dataCorr(groupID,behavIdx);
        [f,x] = ecdf(currData);
        plot(x,f,'LineWidth',1.5,'color',allColors(iV,:))
        
        sigCellsA = data_InModel(groupID,behavIdx);
        for iX = iV:nVar
            sigCellsB = data_InModel(groupID,ismember(Params.BehavVar,varOfInterest{iX}));
            nSig(iX,iV) = sum(sigCellsA & sigCellsB);
        end
        
         [pResults(iV),h,statsResults{iV}] = signrank(currData,shuffCorr(groupID,behavIdx));
    end
    for iV = 1:nVar
        behavIdx = ismember(Params.BehavVar,varOfInterest{iV});
        currData = shuffCorr(groupID,behavIdx);
        [f,x] = ecdf(currData);
        plot(x,f,'k--','LineWidth',.5)
    end
    
    axis([0 1 0 1])
    plotHorizLine(.5,{':r'})
    plotVertLine(nanmedian(reshape(shuffCorr(groupID,cellfun(@(x) find(ismember(Params.BehavVar,x)),varOfInterest)),[],1)),{':r'})
    legend(cat(2,varOfInterest,'Shuffle'),'location','best')
    xlabel('Abs Correlation')
    ylabel('Cummulative Proportion')
    title(sprintf('Relation of mPFC firing to behavior: %s',anatGroup{iA}))
    
    figFN = [figOutputDir,currSubDir,'FullMazeCorrealtions_',anatGroup{iA}];
    fh_data.Position = [100 100 810 620];
    saveas(fh_data,figFN,'svg'); close(fh_data)
    
    
    propSigCells = nSig/sum(groupID);
    fh_data = figure;
    s = imagesc(propSigCells);
    set(s,'alphadata',~isnan(propSigCells));
    caxis([0 1])
    xticks(1:nVar); xticklabels(varOfInterest)
    xtickangle(45)
    yticks(1:nVar); yticklabels(varOfInterest)
    ytickangle(45)
    box off
    colorbar
    title(anatGroup{iA})
    
    figFN = [figOutputDir,currSubDir,'PropCoCorrelations_',anatGroup{iA}];
    fh_data.Position = [100 100 810 620];
    saveas(fh_data,figFN,'svg'); close(fh_data)
    
    
    validCond = cellfun(@(x) find(ismember(Params.BehavVar,x)),varOfInterest);
    fh_data = figure;
    histogram(sum(data_InModel(groupID,validCond),2),-.5:nVar+.5,'normalization','probability');
    box off
    xlabel('Number of sig Variables')
    ylabel('Proportion of Cells')
    title(anatGroup{iA})
    
    figFN = [figOutputDir,currSubDir,'NumberOfSigVar_',anatGroup{iA}];
    saveas(fh_data,figFN,'svg'); close(fh_data)
    
end

%% Plot the correlation of behavioral variables

% Get any behavioral Task Data that may be needed

[FullLapData, FullSessionData, FullIdPhi] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{cat(1,Params.SSN{:})});
[nSess,nLaps,nSites] = size(FullLapData.Decision);
nonLapVar = {'MazeChunk'};

behavData = cell(nVar,1);
for iV = 1:nVar
    if ismember(varOfInterest{iV},nonLapVar)
        behavData{iV} = nan(nSess,nLaps,nSites);
    else
        behavData{iV} = ExtractBehavDataFromFull(FullLapData,FullSessionData,varOfInterest{iV});
    end
end

outputCorr = nan(nVar,nVar,nSess);
for iS = 1:nSess
    for iV = 1:nVar
        for iX = iV+1:nVar
            outputCorr(iX,iV,iS) = nancorr(behavData{iV}(iS,:),behavData{iX}(iS,:));
        end
    end
end

fh_data = figure;
s = imagesc(nanmean(outputCorr,3));
set(s,'alphadata',~isnan(nanmean(outputCorr,3)));
caxis([-1 1])

xticks(1:nVar); xticklabels(varOfInterest)
xtickangle(45)
yticks(1:nVar); yticklabels(varOfInterest)
ytickangle(45)
box off
colorbar
title('Behavior to Behavior Correlation')

figFN = [figOutputDir,currSubDir,'BehavToBehavCorr'];
fh_data.Position = [100 100 810 620];
saveas(fh_data,figFN,'svg'); close(fh_data)

