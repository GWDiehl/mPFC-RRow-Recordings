function mPFC_RevisionAnalysis_LingerZone_eLife(varargin)

figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\LingerZone\';

process_varargin(varargin);


CleanMatlabPlotDefaults
SelectDefaultColorMap(jet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Avg FR in LZ Baseline   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'FiringRate\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

% Avg FR in the linger zone

FNRoot = 'Revision\eLife\PETH\LingerZone\';
dataType = 'General';

GWDColor = loadGWD_RRowColormap;
saveSuffix = '_FR';

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
selectedTime = {[-2 5]};
windowBuffer = [-.15 .15];

dataFN = {'LingerZone_Start'};
selectedEvent = {'LZEntry'};
selectedGroup = {'All'};
plotTitle = {'LingerZone'};

yRange = [-4 1.5];
transitionBin = 0;
binsZ = [-2 0];
binsRes = [.5 3];
peakRate = 6.5;

[subA,subB] = computeSubplotSize(length(anatGroup));

for iX = 1:length(selectedEvent)
    
    fullFN = [FNRoot,dataFN{iX}];
    currEvent = selectedEvent(iX);
    
    % Full slate of groupings
    [fh_data, outputResults, dataCenters, Params] = PlotFullSingleCellResults(fullFN,dataType,currEvent,repmat(selectedGroup,1,length(currEvent)),repmat(plotTitle,1,length(currEvent)),anatGroup,...
        'yRange',yRange,'selectedTime',selectedTime{iX},'windowBuffer',windowBuffer);
    
    zScoreBins = dataCenters >= binsZ(1) & dataCenters <= binsZ(2);
    responseBins = dataCenters>= binsRes(1) & dataCenters <= binsRes(2);
    zData = cell(size(outputResults));
    peakTime = cell(size(outputResults));
    
    fh_PeakTime = figure;   hold on
    fh_FullRaw = figure;
    for iE = 1:numel(outputResults)
        temp = outputResults{iE};
        zData{iE} = (temp - nanmean(temp(:,zScoreBins),2))./nanstd(temp(:,zScoreBins),[],2);
        [tempA,temp] = max(abs(zData{iE}),[],2);
        peakTime{iE} = dataCenters(temp);
        
        figure(fh_PeakTime)
        histogram(peakTime{iE},-2:.25:5,'Normalization','probability','DisplayStyle','stairs','edgecolor',GWDColor.(anatGroup{iE}),'linewidth',1.5)
        plot(nanmean(peakTime{iE}),0,'o','markersize',15,'color',GWDColor.(anatGroup{iE}))
        
        [avgAct,order] = sort(nanmean(outputResults{iE}(:,responseBins),2));
        
        figure(fh_FullRaw)
        subplot(subA,subB,iE)
        s = imagesc(dataCenters,[],outputResults{iE}(order,:));
        set(s,'alphadata',~isnan(outputResults{iE}(order,:)))
        
        plotVertLine(transitionBin,{'r','linewidth',1.5})
        plotHorizLine(find(avgAct>0,1)-.5,{'--r','linewidth',1.5})
        title(anatGroup{iE})
        colorbar
        xlabel('Time (sec)')
        
        if isempty(peakRate)
            temp = prctile(abs(outputResults{iE}(:)),95);
        else
            temp = peakRate;
        end
        caxis([-temp temp])
    end
    figure(fh_PeakTime)
    plotVertLine(0,{'--k'})
    ylabel('Prop Cell')
    xlabel('Peak Response Time (sec)')
    
    nCells = cellfun(@(x) size(x,2),peakTime);
    fullDataList = arrayfun(@(x) repmat(anatGroup(x),1,nCells(x))',1:numel(outputResults),'UniformOutput',0);
    fullDataList = cat(1,fullDataList{:});
    
    fullOutputData = cat(2,peakTime{:})';
    
    [p,ANOVA,ANOVA_stats] = anova1(fullOutputData,fullDataList,'off');
    multiCompare = multcompare(ANOVA_stats,'display','off');
    
    dmPFC = cellfun(@(x) ismember(x,{'Acc_TE' 'dPL_TE'}),fullDataList);
    vmPFC = cellfun(@(x) ismember(x,{'vPL_TE' 'IL_TE'}),fullDataList);
    
    [p,~,stats] = ranksum(fullOutputData(dmPFC),fullOutputData(vmPFC));
    
    figFN = [figOutputDir,currSubDir,sprintf('%s_%s_MainData',selectedEvent{iX},plotTitle{1}),saveSuffix];
    saveas(fh_data,figFN,'svg'); close(fh_data)
    
    figFN = [figOutputDir,currSubDir,sprintf('%s_%s_PeakResponse',selectedEvent{iX},plotTitle{1}),saveSuffix];
    fh_PeakTime.Position = [100 100 1150 620];
    saveas(fh_PeakTime,figFN,'svg'); close(fh_PeakTime)
    
    figFN = [figOutputDir,currSubDir,sprintf('%s_%s_FullFR_RawSorting',selectedEvent{iX},plotTitle{1}),saveSuffix];
    fh_FullRaw.Position = [100 100 1150 620];
    saveas(fh_FullRaw,figFN,'svg'); close(fh_FullRaw)
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Avg FR in LZ By Linger Time   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'FiringRate\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

% Linger Zone FR split by time spent lingering

FNRoot = 'Revision\eLife\PETH\LingerZone\';
dataType = 'General';

GWDColor = loadGWD_RRowColormap;
saveSuffix = '_FR';

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
selectedTime = 8:33;
windowBuffer = [-.5 .5];

dataFN = 'LingerZone_Warped_SplitByTime';
currEvent = {'LZEntry_PrctLingerTime0_5','LZEntry_PrctLingerTime5_10','LZEntry_PrctLingerTime10_15','LZEntry_PrctLingerTime15_20','LZEntry_PrctLingerTime20_25','LZEntry_PrctLingerTime25_30'};
selectedGroup = repmat({'All'},1,length(currEvent));
plotTitle = repmat({'LingerZone'},1,length(currEvent));

yRange = [];
transitionBin = [10 30];
[subA,subB] = computeSubplotSize(length(anatGroup));
allColors = jet(length(currEvent));


fullFN = [FNRoot,dataFN];

% Full slate of groupings
[fh_data, outputResults, dataCenters, Params] = PlotFullSingleCellResults(fullFN,dataType,currEvent,selectedGroup,plotTitle,anatGroup,...
    'yRange',yRange,'selectedTime',selectedTime,'windowBuffer',windowBuffer);
close(fh_data)

fh_data = figure;
for iA = 1:length(anatGroup)
    subplot(subA,subB,iA); hold on
    
    for iE = 1:length(currEvent)
        PlotSingleCellAvg(outputResults{iE,iA},dataCenters,'plotColor',allColors(iE,:))
    end
    title(anatGroup{iA})
    xlabel('Norm Time')
    ylabel('Norm FR')
    
    if ~isempty(yRange)
        ylim(yRange)
    end
    
    xlim([dataCenters(1)+windowBuffer(1) dataCenters(end)+windowBuffer(2)])
    plotVertLine(transitionBin,{'--k'})
    plotHorizLine(0,{'--k'})
end

figFN = [figOutputDir,currSubDir,'LZEntry_LingerTime_MainData',saveSuffix];
fh_data.Position = [100 100 1150 620];
saveas(fh_data,figFN,'svg'); close(fh_data)


%%

GWDColor = loadGWD_RRowColormap;
quantWindows = [1,10;10,20;20,30;30,40];
quantNames = {'PreRZ' 'EarlyRZ' 'LateRZ' 'PostRZ'};
eventData = [2.5:5:27.5];


cellIDs = cat(1,Params.CellID{:});
groupID = groupPFCCellsByDepth(cellIDs,anatGroup);
ratID = cellfun(@(x) x(1:4),cellIDs,'UniformOutput',0);
[allRats,~,ratIdx] = unique(ratID,'stable');
nRats = max(ratIdx);
nAnat = length(anatGroup);

nQuants = [0 size(quantWindows,1)];
totalQuants = sum(nQuants);

nEvents = length(currEvent);
meanResponse = nan(nEvents,nAnat,nRats,totalQuants);
netCorr = nan(nAnat,totalQuants);
netError = nan(nAnat,totalQuants,2);
netSig = nan(nAnat,totalQuants);

fh_data = figure; hold on
for iA = 1:nAnat
    for iQ = 1:nQuants(iX+1)
        selectedBins = dataCenters >=quantWindows(iQ,1) & dataCenters <= quantWindows(iQ,2);
        
        iQQ = sum(nQuants(1:iX))+iQ;
        for iR = 1:nRats
            selectedCells = ratIdx(groupID(:,iA)) == iR;
            for iE = 1:nEvents
                tempResponse = nanmean(outputResults{iE,iA}(selectedCells,:),1);
                meanResponse(iE,iA,iR,iQQ) = nanmean(tempResponse(selectedBins));
            end
        end
        tempA = reshape(repmat(eventData,1,nRats),[],1);
        tempB = reshape(meanResponse(:,iA,:,iQQ),[],1);
        valid = ~isnan(tempA+tempB);
        
        [~,~,RL,RU] = corrcoef(tempA(valid),tempB(valid));
        netError(iA,iQQ,1) = RL(2); netError(iA,iQQ,2) = RU(2);
        
        [netCorr(iA,iQQ),netSig(iA,iQQ)] = nancorr(tempA(valid),tempB(valid));
    end
    
%     ShadedErrorbar(1:totalQuants, netCorr(iA,:), [],
%     'L',-netError(iA,:,1)+netCorr(iA,:),'U',netError(iA,:,2)-netCorr(iA,:),'color',GWDColor.(anatGroup{iA}),'lineWidth',1,'marker','none');
    plot(1:totalQuants,netCorr(iA,:),'--o','color',GWDColor.(anatGroup{iA}));
end

ylim([-.5 .5])
xlim([.5 totalQuants+.5])
xticks(1:totalQuants);xticklabels(quantNames)
plotHorizLine(0,{'--k'})
ylabel('Correlation')

figFN = [figOutputDir,currSubDir,sprintf('FRtoLingerTime_Corr')];
saveas(fh_data,figFN,'svg'); close(fh_data)

