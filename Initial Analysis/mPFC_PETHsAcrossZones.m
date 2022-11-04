function mPFC_PETHSAcrossZones

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
outputRoot = 'Revision\eLife\PETH\';

analysisType = 'RawFR';
dt = .1 ;
nShuffles = 30;


%% Offer Zone FR

subDir = 'OfferZone';
EventLimits = {'','Skip','Accept'};

outputFN = {'OfferZone_Warped'};
RRowEvent = {'OZEntry'};

timingType = 'WarpedTime';

timeWindow = [10 20 10];
eventTimeBounds = [-1 0 1 2];
boundOverride = struct;
boundOverride.Var = 'Accept';
boundOverride.Val = [.5 1.5];
boundOverride.Ref = [nan nan nan 3];
boundOverride.Dur = [nan nan nan 1];

for iX = 1:length(outputFN)
    
    [Results, Params] = CreatePETH_MI_Baseline('RRowEvent',RRowEvent(iX),'EventLimits',EventLimits,...
        'timingType',timingType,'analysisType',analysisType,'eventTimeBounds',eventTimeBounds,'boundOverride',boundOverride,...
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

fprintf('Done with Offer Zone \n')

%% Wait Zone FR Split by Delay and Value

subDir = 'WaitZone';
EventLimits = {{'' 'Delay1,5' 'Delay5,10' 'Delay10,15' 'Delay15,20' 'Delay20,25' 'Delay25,30' 'Delay1,3' 'Delay4,6' 'Delay7,9' 'Delay10,12' 'Delay13,15' 'Delay16,18' 'Delay19,21' 'Delay22,24'}...
    {'' 'Value-15,-10' 'Value-10,-5' 'Value-5,0' 'Value0,5' 'Value5,10' 'Value10,15'}};
limitSuffix = {'_SplitByDelay' '_SplitByValue'};

outputFN = {'WaitZone_Start' 'WaitZone_End'};
RRowEvent = {'WZEntry' 'WZExit'};

timingType = 'RawTime';

timeWindow = {[-2 10] [-10 2]};
eventTimeBounds = [-1 1];
boundOverride = [];

for iX = 1:length(outputFN)
    for iY = 1:length(limitSuffix)
        
        [Results, Params] = CreatePETH_MI_Baseline('RRowEvent',RRowEvent(iX),'EventLimits',EventLimits{iY},...
            'timingType',timingType,'analysisType',analysisType,'eventTimeBounds',eventTimeBounds,'boundOverride',boundOverride,...
            'timeWindow',timeWindow{iX},'dt',dt,'nShuffles',nShuffles);
        
        fullOutputDir = fullfile(RootDir,AnalysisDir,outputRoot,subDir);
        if ~exist(fullOutputDir,'dir')
            mkdir(fullOutputDir)
        end
        tempOutput = struct;
        tempOutput.Results = Results;
        tempOutput.Params = Params;
        save([fullOutputDir,'\',outputFN{iX},limitSuffix{iY}],'-struct','tempOutput')
    end
end

fprintf('Done with Wait Zone \n')

%% Linger Zone FR at the start

subDir = 'LingerZone';
EventLimits = {''};

outputFN = {'LignerZone_Start'};
RRowEvent = {'LZEntry'};

timingType = 'RawTime';

timeWindow = {[-2 10]};
eventTimeBounds = [-1 1];
boundOverride = [];

for iX = 1:length(outputFN)
        
        [Results, Params] = CreatePETH_MI_Baseline('RRowEvent',RRowEvent(iX),'EventLimits',EventLimits,...
            'timingType',timingType,'analysisType',analysisType,'eventTimeBounds',eventTimeBounds,'boundOverride',boundOverride,...
            'timeWindow',timeWindow{iX},'dt',dt,'nShuffles',nShuffles);
        
        fullOutputDir = fullfile(RootDir,AnalysisDir,outputRoot,subDir);
        if ~exist(fullOutputDir,'dir')
            mkdir(fullOutputDir)
        end
        tempOutput = struct;
        tempOutput.Results = Results;
        tempOutput.Params = Params;
        save([fullOutputDir,'\',outputFN{iX}],'-struct','tempOutput')
end

fprintf('Done with Linger Zone \n')

%% Linger Zone warped and split by linger time

subDir = 'LingerZone';
EventLimits = {'' 'LingerTime3,5' 'LingerTime5,7' 'LingerTime7,9' 'LingerTime9,11'...
    'PrctLingerTime0,5' 'PrctLingerTime5,10' 'PrctLingerTime10,15' 'PrctLingerTime15,20' 'PrctLingerTime20,25' 'PrctLingerTime25,30'};

outputFN = {'LingerZone_Warped_SplitByTime'};
RRowEvent = {'LZEntry'};

timingType = 'WarpedTime';

timeWindow = [10 20 10];
eventTimeBounds = [-1 0 1 2];
boundOverride = struct;
boundOverride.Var = 'Earn';
boundOverride.Val = [.5 1.5];
boundOverride.Ref = [2 nan nan nan];
boundOverride.Dur = [-1 nan nan nan];

for iX = 1:length(outputFN)
    
    [Results, Params] = CreatePETH_MI_Baseline('RRowEvent',RRowEvent(iX),'EventLimits',EventLimits,...
        'timingType',timingType,'analysisType',analysisType,'eventTimeBounds',eventTimeBounds,'boundOverride',boundOverride,...
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

fprintf('Done with Linger Zone: Linger TIme \n')


