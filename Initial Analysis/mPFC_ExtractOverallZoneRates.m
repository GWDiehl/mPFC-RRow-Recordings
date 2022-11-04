function mPFC_ExtractOverallZoneRates

% Shell function to compute the avg firing rate in various zones of the
% task. Compute for each cell along with shuffled data for each cell to
% normalize against.

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
outputRoot = 'Revision\eLife\AvgEventRate\';

nShuffles = 30;


fullOutputDir = fullfile(RootDir,AnalysisDir,outputRoot);
if ~exist(fullOutputDir,'dir')
    mkdir(fullOutputDir)
end


%% Look at the avg firing rate throughout a period of the behavioral task


eventStart = {'LapStart' 'OZEntry' 'WZEntry' 'LZEntry' 'AZEntry'};
eventEnd = {'LapEnd' 'OZExit' 'WZExit' 'LZExit' 'AZExit'};
EventLimits = {'' '' '' '' ''};

[Results,Params] = ExtractAvgRateBetweenEvents('eventStart',eventStart,'eventEnd',eventEnd,'EventLimits',EventLimits,'nShuffles',nShuffles);

outputFN = 'GeneralEvents_FullWindow';
save([fullOutputDir,'\',outputFN],'Results','Params')

