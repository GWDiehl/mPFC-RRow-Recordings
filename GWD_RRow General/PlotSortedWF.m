
function [fh1,fh2,sortedWF,peakChan] = PlotSortedWF(S_FN,channelMap_FN,varargin)



samplingRate = 30; % Sampling rate in kHz

plotRange = 3;
baselineColor = [.3 .3 .3 .3];

TFileRoot = '.t'; % What is the suffix on the spike time files
WaveFormSuffix = '_WF';

[RootDir, ~, ~, ~, ~, ~, SpikesDir] = ID_GWD_Dirs;

plottingType = 'Local';


process_varargin(varargin);

%%

SSN = S_FN(1:15);
RatID = S_FN(1:4);
WF_FN = [S_FN TFileRoot(2:end) WaveFormSuffix];

fullWF_FN = fullfile(RootDir,SpikesDir,RatID,SSN,WF_FN);

[sortedWF,peakChan] = LoadAndSortProbeWF({fullWF_FN},{channelMap_FN});
sortedWF = sortedWF{1};
[nSamples,nChan] = size(sortedWF);


chanWindow = peakChan-plotRange:peakChan+plotRange;
validWindow = chanWindow>=1 & chanWindow<=nChan;
idxOfInterest = false(nChan,1);
idxOfInterest(chanWindow(validWindow)) = 1;

switch plottingType
    case 'Local'
        colorSet = jet(plotRange*2+1);
        colorSet = colorSet(validWindow,:);
        colorSet(:,4) = 1;
        plotColors = repmat(baselineColor,nChan,1);
        plotColors(idxOfInterest,:) = colorSet;
    case 'Global'
        plotColors = jet(nChan);
        plotColors(:,4) = 1;
        plotColors(~idxOfInterest,4) = baselineColor(4);
end


fh1 = figure;
hold on
for iX = 1:nChan
    
    plot([1:nSamples]/samplingRate,sortedWF(:,iX),'color',plotColors(iX,:))
end
plot([1:nSamples]/samplingRate,sortedWF(:,peakChan),'--k','linewidth',1)
xlabel('Time (ms)')
ylabel('Voltage (uV)')
title(S_FN)

fh2 = figure;
imagesc([1:nSamples]/samplingRate,1:nChan,sortedWF')
xlabel('Time (ms)')
ylabel('Channel Number: Linear Position')
colorbar
title(S_FN)

