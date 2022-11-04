function mPFC_LinearRegressionToBehavior


%% Compute the linear regression between cell firing rate and behavioral variables

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
outputRoot = 'Revision\eLife\SingleCellRegression\';

% Make sure the directory exists
if ~exist([RootDir,AnalysisDir,outputRoot],'dir')
    mkdir([RootDir,AnalysisDir,outputRoot])
end


analysisType = 'ByChunk';
nShuffles = 30;

if isempty(gcp('nocreate'))
    parpool(6)
end
pctRunOnAll warning off

outputFN = 'BasicRelationToBehavior';

analysisWindow = {'OZEnter' 'WZEnter' 'LZEnter' 'AZEnter'; 'OZExit' 'WZExit' 'LZExit' 'AZEnter'};
behavVar = {'SessionTime' 'MazeChunk' 'Skip' 'Value' 'Delay' 'Rank' 'OZRunSpeed' 'Quit' 'ReactionTime' 'LingerTime' 'EconomicChoice' 'RandomLap'};

[Results, Params] = RRow_RegressionToBehavior('analysisType',analysisType,...
    'analysisWindow',analysisWindow,'behavVar',behavVar,'analysesToRun',{'Correlation' 'MutualInfo'},'nShuffles',nShuffles);

save([RootDir,AnalysisDir,outputRoot,outputFN],'Results','Params')

fprintf('\n Done with %s \n\n','FullSession Data')

