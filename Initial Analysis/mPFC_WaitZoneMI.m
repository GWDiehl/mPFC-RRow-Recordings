function mPFC_WaitZoneMI


%% Compute the mutual information between cell spiking and behaviors (Wait Zone)

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
outputRoot = 'Revision\eLife\MutualInfo\';

analysisType = 'MutualInfo';
dt = .1 ;
nShuffles = 30;

%% Wait Zone Quit/Earn and the Start/End

subDir = 'WaitZone';

outputFN = {'WaitZone_Start' 'WaitZone_End'};
RRowEvent = {'WZEntry' 'WZExit'};

EventLimits = {''};
MI_Params = struct;
MI_Params.Var = 'Quit';
MI_Params.BehavBins = [nan nan 2]; % Start/End/nBins
MI_Params.SpkBins = [1 99 5]; % StartPrctile/EndPrctile/nBins
MI_Params.UseLapData = 1;

timingType = 'RawTime';
timeWindow = [-2 10];
eventTimeBounds = [-1 1];
boundOverride = [];

for iX = 2:length(outputFN)
    
    [Results, Params] = CreatePETH_MI_Baseline('RRowEvent',RRowEvent(iX),'EventLimits',EventLimits,...
        'timingType',timingType,'analysisType',analysisType,'MI_Params',MI_Params,...
        'eventTimeBounds',eventTimeBounds,'boundOverride',boundOverride,...
        'timeWindow',timeWindow,'dt',dt,'nShuffles',nShuffles);
    
    fullOutputDir = fullfile(RootDir,AnalysisDir,outputRoot,subDir);
    if ~exist(fullOutputDir,'dir')
        mkdir(fullOutputDir)
    end
    tempOutput = struct;
    tempOutput.Results = Results;
    tempOutput.Params = Params;
    save([fullOutputDir,'\',outputFN{iX}],'-struct','tempOutput')    
end



