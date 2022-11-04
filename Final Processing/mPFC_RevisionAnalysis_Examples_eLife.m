function mPFC_RevisionAnalysis_Examples_eLife(varargin)

figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\';

process_varargin(varargin);

CleanMatlabPlotDefaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Specific  examples of firing to behavior   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'ExampleCells\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

%%


behavVar = {'Delay' 'SessionTime' 'Rank' 'Value'}; %
behavDisc = {[.5 30.5 30] [0 3600 30] [.5 4.5 4] [-15 15 30]};

selectedCells = {'R530-2019-06-19-Si02_02';'R537-2019-12-11-Si01_04';'R535-2019-07-25-Si02_27';'R535-2019-07-19-Si02_13';'R542-2019-07-25-Si02_32'};

xRange = [7 33];
timeWindowLimit = [8:33];

% Load in the stepwise regression data and get cellIDs
[RootDir, ~, BehavDir, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;

eventsToUse = 1;
dataRootName = 'Revision\eLife\AvgEventRate\';
dataFN = 'GeneralEvents_FullWindow';


nBehav = length(behavVar);
nCells = length(selectedCells);

temp = load([RootDir,AnalysisDir,dataRootName,dataFN]);
Params = temp.Params;
Results = temp.Results;

FRData = Results.LapStart_LapEnd.AvgResults;


% Get any behavioral Task Data that may be needed
[FullLapData, FullSessionData, FullIdPhi] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{cat(1,Params.SSN{:})});

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
        behavDisc{iB} = [.5 nEvents+.5 nEvents];
    end
end

%%
[subA,subB] = computeSubplotSize(nBehav);

for iC = 1:nCells
    currSubDir = ['ExampleCells\',selectedCells{iC},'\'];
    
    if ~exist([figOutputDir,currSubDir],'dir')
        mkdir([figOutputDir,currSubDir]);
    end
    
    sessIdx = find(cellfun(@(x) ismember(selectedCells{iC},x),Params.CellID));
    cellIdx = find(ismember(Params.CellID{sessIdx},selectedCells{iC}));
    
    FRProfile = nan(1,0);
    behavProfile = nan(nBehav,0);
    behavProfile_Raw = nan(nBehav,0);
    
    nSamples = size(FRData{sessIdx},2);
    lapIdx = Results.LapStart_LapEnd.LapIdx{sessIdx};
    FRProfile = cat(2,FRProfile,FRData{sessIdx}(cellIdx,:));
    tempBehav = nan(nBehav,nSamples,2);
    for iB = 1:nBehav
        if strcmp(behavVar{iB},'MazeChunk')
            tempBehav(iB,:) = iE;
        else
            temp = reshape(binnedBehav{iB}(sessIdx,:),nLaps,nSites);
            tempBehav(iB,:,1) = arrayfun(@(x) temp(behavIdxKey(sessIdx,:)==x),lapIdx);
            
            temp = reshape(behavData{iB}(sessIdx,:),nLaps,nSites);
            tempBehav(iB,:,2) = arrayfun(@(x) temp(behavIdxKey(sessIdx,:)==x),lapIdx);
        end
    end
    behavProfile = cat(2,behavProfile,tempBehav(:,:,1));
    behavProfile_Raw = cat(2,behavProfile_Raw,tempBehav(:,:,2));
    
    fh_data = figure;
    for iB = 1:nBehav
        subplot(subA,subB,iB); hold on
        plot(behavProfile_Raw(iB,:)+(rand(1,nSamples)-.5)*.5,FRProfile,'r.')
        
        temp = num2cell(behavDisc{iB}); temp{end} = temp{end}+1;
        meanData = arrayfun(@(x) nanmean(FRProfile(:,behavProfile(iB,:)==x),2),1:behavDisc{iB}(3));
        errorData = arrayfun(@(x) nanstderr(FRProfile(:,behavProfile(iB,:)==x)),1:behavDisc{iB}(3));
        ShadedErrorbar(computeBinCenters(linspace(temp{:})),meanData,errorData,'lineWidth',1,'marker','none');
        errorbar(computeBinCenters(linspace(temp{:})),meanData,errorData,'b.');
        
        title(sprintf('%s: Corr of %0.2f',behavVar{iB},nancorr(FRProfile,behavProfile_Raw(iB,:))))
        xlim([temp{1},temp{2}])
        xlabel(behavVar{iB})
        ylabel('Avg Firing Rate')
    end
    figFN = [figOutputDir,currSubDir,'RelatingToBehav'];
    fh_data.Position = [100 100 810 620];
    saveas(fh_data,figFN,'svg'); close(fh_data)
    
    if ismember('Rank',behavVar)
        fh_data = figure;
        varIdx = ismember(behavVar,'Rank');
        meanData = arrayfun(@(x) nanmean(FRProfile(:,behavProfile(varIdx,:)==x),2),1:behavDisc{varIdx}(3));
        imagesc(meanData)
        xlabel('Rank')
        temp = reshape(binnedBehav{varIdx}(sessIdx,:),nLaps,nSites);
        temp = temp(1,:);
        [~,order] = sort(temp,2);
        title(sprintf('Ranked Site Num: %d, %d, %d, %d',order(1),order(2),order(3),order(4)))
        colorbar
        
        figFN = [figOutputDir,currSubDir,'RelatingToBehav_RankHeatmap'];
        saveas(fh_data,figFN,'svg'); close(fh_data)
    end    
end

rawFN = [RootDir,AnalysisDir,'\Revision\eLife\PETH\OfferZone\OfferZone_Warped'];
rawData = load(rawFN);

events = {'OZEntry' 'OZEntry_Accept' 'OZEntry_Skip'};
fh_data = PlotCellExampleResponse(selectedCells,events,rawData,'General','timeWindowLimit',timeWindowLimit,'xRange',xRange);
figFN = [figOutputDir,'ExampleCells\PETH_OZ'];
fh_data.Position = [100 100 810 620];
saveas(fh_data,figFN,'svg'); close(fh_data)

