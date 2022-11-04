function mPFC_RevisionAnalysis_Behavior_eLife(varargin)


figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\';


process_varargin(varargin);


CleanMatlabPlotDefaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Behavioral Data   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'Behavior\';
if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

[rootDir, DataDefName,BehavDir] = ID_GWD_Dirs;

dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

SSNs = dataDef.SSN(OI);

% Get any behavioral Task Data that may be needed
[FullLapData, FullSessionData, FullIdPhi] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});


%% Threhsold example

SSN_OI =  'R536-2019-09-18';

sessionIdx = ismember(SSN_OI,dataDef.SSN(OI));

waitZones = find(ismember(FullLapData.ZoneID(sessionIdx,1,:),'WaitZone'));
offerZones = find(ismember(FullLapData.ZoneID(sessionIdx,1,:),'OfferZone'));

delay = squeeze(FullLapData.ZoneDelay(sessionIdx,:,offerZones));
accept = squeeze(FullLapData.AcceptOffer(sessionIdx,:,offerZones)==1);
quit = squeeze(FullLapData.QuitOffer(sessionIdx,:,waitZones)==1);

fh = plotRRowChoices(delay,accept,quit,'thresholdType','Heaviside');

figFN = [figOutputDir,currSubDir,'ThresholdExample'];
saveas(fh,figFN,'svg'); close(fh)


%% Overall proportion of Accept/Skip/Quit/Earned Laps

%%% Note: This does not account for conditional relationships/probabiltiies
%%% and thus will sum to >1

allRats = unique(dataDef.RatID(OI),'stable');
nGrps = length(allRats);
choices = {'Skip' 'Accept' 'Quit_Global' 'Rewarded' 'Quit' 'Earn'};
ChoiceData = nan(nGrps,length(choices),2);

fh = figure; hold on
grpColors = jet(nGrps);
for iC = 1:length(choices)
    [fh_temp,~,~,ChoiceData(:,iC,:)] = Plot2DBehavData_RRow({'Null' choices{iC}},'groupingFunc',@nansum,'normalizeResults','LapCount');
    close(fh_temp);
    errorbar(iC,nanmean(ChoiceData(:,iC,2)),nanstd(ChoiceData(:,iC,2))/sqrt(nGrps),'k','linewidth',1)
    plot(iC,nanmean(ChoiceData(:,iC,2)),'xk','linewidth',1)
end
figure(fh)
for iG = 1:nGrps
    plot(ChoiceData(iG,:,2),'o','color',grpColors(iG,:));
end
axis([.5 length(choices)+.5 -.05 1])
ylabel('Proportion of All laps')
tempLabel = cellfun(@(x) [x,'    '],choices,'UniformOutput',0);
xlabel([tempLabel{:}])
figFN = [figOutputDir,currSubDir,'OverallChoices'];
saveas(fh,figFN,'svg'); close(fh)


%% Accept/Quit Probabilities as a function of offer
% Also RT bc uses the same code

% Accept Prob
fh = Plot2DBehavData_RRow({'Accept' 'Value'});
figure(fh); 
plotVertLine(0,{'--k'})
legend(unique(dataDef.RatID(OI),'stable'),'location','best');

figFN = [figOutputDir,currSubDir,'AcceptProb'];
saveas(fh,figFN,'svg'); close(fh)

% Quit Prob
fh = Plot2DBehavData_RRow({'Quit' 'Value'});
figure(fh);plotVertLine(0,{'--k'})
legend(unique(dataDef.RatID(OI),'stable'),'location','best');

figFN = [figOutputDir,currSubDir,'QuitProb'];
saveas(fh,figFN,'svg');  close(fh)

% Reaction time
fh = Plot2DBehavData_RRow({'ReactionTime' 'Value'});
figure(fh);plotVertLine(0,{'--k'})
legend(unique(dataDef.RatID(OI),'stable'),'location','best');

figFN = [figOutputDir,currSubDir,'ReactionTime'];
saveas(fh,figFN,'svg');  close(fh)

% Quit Time
fh = Plot2DBehavData_RRow({'QuitTime' 'Delay'});
plotIDLine({':r'})
figFN = [figOutputDir,currSubDir,'QuitTime_ByDelay'];
saveas(fh,figFN,'svg');  close(fh)
[fh,~,~,quitData,binEdges] = Plot2DBehavData_RRow({'QuitTime' 'Value'});
figFN = [figOutputDir,currSubDir,'QuitTime_ByValue'];
saveas(fh,figFN,'svg');  close(fh)
valueBins = computeBinCenters(binEdges{2});
[QT_Corr, QT_pVal] = nancorr(repmat(valueBins,size(quitData,1),1),quitData);

quitTimes = ExtractBehavDataFromFull(FullLapData,[],'QuitTime');
fh = figure; histogram(quitTimes(:),0:.75:27);
box off
xlabel('Quit Time (sec)')
ylabel('Number of Quits')
figFN = [figOutputDir,currSubDir,'QuitTime'];
saveas(fh,figFN,'svg');  close(fh)

% Linger Time
fh = Plot2DBehavData_RRow({'LingerTime' 'Value'});
figFN = [figOutputDir,currSubDir,'LingerTime_ByValue'];
saveas(fh,figFN,'svg');  close(fh)
fh = Plot2DBehavData_RRow({'LingerTime' 'Delay'});
figFN = [figOutputDir,currSubDir,'LingerTime_ByDelay'];
saveas(fh,figFN,'svg');  close(fh)

[lingerTime,baseEdges] = ExtractBehavDataFromFull(FullLapData,[],'LingerTime');
fh = figure; histogram(lingerTime(:),2.5:.5:35);
box off
xlabel('Linger Time (sec)')
ylabel('Number of Earned Rewards')
figFN = [figOutputDir,currSubDir,'LingerTime'];
saveas(fh,figFN,'svg');  close(fh)

% Reaction time by Choice
[fh,~,~,acceptData] = Plot2DBehavData_RRow({'AcceptRT' 'Value'},'plotGrps',0);
[fh1,~,~,skipData,binEdges] = Plot2DBehavData_RRow({'SkipRT' 'Value'},'plotGrps',0);
valueBins = computeBinCenters(binEdges{2});

[p, ~, stats] = signrank(acceptData(:),skipData(:));


close(fh1); figure(fh);hold on
dataMean = nanmean(skipData,1);
dataSEM = nanstd(skipData,[],1)./sqrt(sum(~isnan(skipData),1));
errorbar(valueBins,dataMean,dataSEM,'r','LineWidth',1.5)
box off
plotVertLine(0,{'--k'})
legend({'Accept' 'Skip'})
ylabel('Reaction Time')

figFN = [figOutputDir,currSubDir,'ReactionTime_Choice'];
saveas(fh,figFN,'svg');  close(fh)

[RT_Corr, RT_pVal] = nancorr(repmat(valueBins,size(acceptData,1),1),acceptData);
posVal = valueBins>0;
[RT_Corr, RT_pVal] = nancorr(repmat(valueBins(posVal),size(acceptData,1),1),acceptData(:,posVal));
negVal = valueBins<0;
[RT_Corr, RT_pVal] = nancorr(repmat(valueBins(negVal),size(acceptData,1),1),acceptData(:,negVal));

acceptRT = ExtractBehavDataFromFull(FullLapData,[],'AcceptRT');
acceptRT = acceptRT(~isnan(acceptRT));
skipRT = ExtractBehavDataFromFull(FullLapData,[],'SkipRT');
skipRT = skipRT(~isnan(skipRT));

[p,~,stats] = ranksum(acceptRT,skipRT);


%% Threshold quantifications indicating stable and subjective

allRats = unique(FullSessionData.RatID,'stable');
groupingID = cellfun(@(x) find(ismember(allRats,x)),FullSessionData.RatID);
thresholds = ExtractBehavDataFromFull(FullLapData,FullSessionData,'Threshold');

[avgThresh, threshSpread, threshCorr] = PlotThresholdMetrics(reshape(thresholds(:,1,:),size(thresholds,1),size(thresholds,3)),groupingID);

figFN = [figOutputDir,currSubDir,'AvgThrehsold'];
saveas(avgThresh,figFN,'svg'); close(avgThresh)

figFN = [figOutputDir,currSubDir,'ThresholdVar'];
saveas(threshSpread,figFN,'svg'); close(threshSpread)

figFN = [figOutputDir,currSubDir,'ThresholdCorr'];
saveas(threshCorr,figFN,'svg'); close(threshCorr)



