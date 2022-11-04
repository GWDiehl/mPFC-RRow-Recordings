function mPFC_RevisionAnalysis_WaitZone_eLife(varargin)


figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\WaitZone\';
process_varargin(varargin);

CleanMatlabPlotDefaults


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Mutual Information for single cells   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'MutualInfo\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

% MI for Quit Choice

FNRoot = 'Revision\eLife\MutualInfo\WaitZone\';
dataType = 'General';

saveSuffix = '_MI';

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
decisionWindow = {[0 2] [-2 0]};;
selectedTime = {[-2 5] [-5 2]};
windowBuffer = [-.15 .15];

dataFN = {'WaitZone_Start' 'WaitZone_End'};
selectedEvent = {'WZEntry' 'WZExit'};
selectedGroup = {'All'};
plotTitle = {'Choice'};

yRange = [-.01 .03];
sigYPos = linspace(.02,.025,length(anatGroup))';

for iX = 1:length(selectedEvent)
    
    fullFN = [FNRoot,dataFN{iX},'\'];
    
    % Full slate of groupings
    [fh_data, outputResults, dataCenters, Params] = PlotFullSingleCellResults(fullFN,dataType,selectedEvent(iX),selectedGroup,plotTitle,anatGroup,...
        'yRange',yRange,'selectedTime',selectedTime{iX},'windowBuffer',windowBuffer);
    ylabel('Norm. Mutual Info')
    
    decisionBins = dataCenters >=decisionWindow{iX}(1) & dataCenters <= decisionWindow{iX}(2);
    % Stats w/ significance for the full data
    [rmANOVA,pVals,sigTime,ANOVA,ANOVA_stats,multiCompare] = ComputeDecisionStats_mPFC(outputResults,anatGroup,decisionBins);
    
    % Plot the Sig Fig
    [fh_sig] = PlotSingleCellSignificance(sigTime,Params,dataCenters,anatGroup,'xRange',[dataCenters(1)+windowBuffer(1) dataCenters(end)+windowBuffer(2)]);
    title(plotTitle{1})
    
     % Overlay on the original Avg data
        [fh_data] = PlotSingleCellSignificance_OnAvg(fh_data,sigTime,dataCenters,anatGroup,sigYPos);
    
    set(fh_data,'name',sprintf('MI For %s',plotTitle{1}));
    set(fh_sig,'name',sprintf('MI For %s',plotTitle{1}));
    
    figFN = [figOutputDir,currSubDir,sprintf('%s_%s_MainData',selectedEvent{iX},plotTitle{1}),saveSuffix];
    saveas(fh_data,figFN,'svg'); close(fh_data)
    figFN = [figOutputDir,currSubDir,sprintf('%s_SigVal',selectedEvent{iX},plotTitle{1}),saveSuffix];
    saveas(fh_sig,figFN,'svg'); close(fh_sig)
end

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Avg FR in WZ Split by Value/Delay   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currSubDir = 'FiringRate\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

% Avg FR split by Value

FNRoot = 'Revision\eLife\PETH\WaitZone\';
FN_Suffix = {'_SplitByValue','_SplitByDelay'};
dataType = 'General';

saveSuffix = '_FR';

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
selectedTime = {[-2 5] [-5 2]};
windowBuffer = [-.15 .15];

dataFN = {'WaitZone_Start' 'WaitZone_End'};
selectedEvent = {'WZEntry' 'WZExit'};
eventSplits = {{'Value_15__10','Value_10__5','Value_5_0','Value0_5','Value5_10','Value10_15'}...
    {'Delay1_5','Delay5_10','Delay10_15','Delay15_20','Delay20_25','Delay25_30'}};
selectedGroup = {'All'};
plotTitle = {'ByValue' 'ByDelay'};

yRange = [-3 3.5];
transitionBin = 0;

[subA,subB] = computeSubplotSize(length(anatGroup));

fullOutputData = cell(length(FN_Suffix),length(selectedEvent));
fullOutputCenters = cell(length(FN_Suffix),length(selectedEvent));

for iM = 1:length(FN_Suffix)
    allColors = jet(length(eventSplits{iM}));
    for iX = 1:length(selectedEvent)
        
        fullFN = [FNRoot,dataFN{iX},FN_Suffix{iM}];
        currEvent = cellfun(@(x) [selectedEvent{iX},'_',x],eventSplits{iM},'UniformOutput',0);
        
        % Full slate of groupings
        [fh_data, fullOutputData{iM,iX}, fullOutputCenters{iM,iX}, Params] = PlotFullSingleCellResults(fullFN,dataType,currEvent,repmat(selectedGroup,1,length(currEvent)),repmat(plotTitle,1,length(currEvent)),anatGroup,...
            'yRange',yRange,'selectedTime',selectedTime{iX},'windowBuffer',windowBuffer);
        close(fh_data);
        
        fh_data = figure;
        for iA = 1:length(anatGroup)
            subplot(subA,subB,iA); hold on
            
            for iE = 1:length(currEvent)
                PlotSingleCellAvg(fullOutputData{iM,iX}{iE,iA},fullOutputCenters{iM,iX},'plotColor',allColors(iE,:))
            end
            title(anatGroup{iA})
            xlabel('Time (sec)')
            ylabel('Norm FR')
            
            if ~isempty(yRange)
                ylim(yRange)
            end
            
            xlim([fullOutputCenters{iM,iX}(1)+windowBuffer(1) fullOutputCenters{iM,iX}(end)+windowBuffer(2)])
            plotVertLine(transitionBin,{'--k'})
            plotHorizLine(0,{'--k'})
        end
        
        figFN = [figOutputDir,currSubDir,sprintf('%s_%s_MainData',selectedEvent{iX},plotTitle{iM}),saveSuffix];
        fh_data.Position = [100 100 1150 620];
        saveas(fh_data,figFN,'svg'); close(fh_data)
    end
end

%%

GWDColor = loadGWD_RRowColormap;
quantWindows = {[-1 0;0 3] [-3 0;0 1]};
quantNames = {'PreWZ' 'EarlyWZ' 'LateWZ' 'PostWZ'};
eventData = {[-12.5:5:12.5]' [2.5:5:27.5]'};
yRange = [-.8 .15; -.15 .55];

cellIDs = cat(1,Params.CellID{:});
groupID = groupPFCCellsByDepth(cellIDs,anatGroup);
ratID = cellfun(@(x) x(1:4),cellIDs,'UniformOutput',0);
[allRats,~,ratIdx] = unique(ratID,'stable');
nRats = max(ratIdx);
nAnat = length(anatGroup);

nQuants = [0 cellfun(@(x) size(x,1),quantWindows)];
totalQuants = sum(nQuants);

for iM = 1:length(FN_Suffix)
    nEvents = length(eventSplits{iM});
    meanResponse = nan(nEvents,nAnat,nRats,totalQuants);
    netCorr = nan(nAnat,totalQuants);
    netError = nan(nAnat,totalQuants,2);
    netSig = nan(nAnat,totalQuants);
    
    fh_data = figure; hold on
    for iA = 1:nAnat
        for iX = 1:length(selectedEvent)
            for iQ = 1:nQuants(iX+1)
                selectedBins = fullOutputCenters{iM,iX} >=quantWindows{iX}(iQ,1) & fullOutputCenters{iM,iX} <= quantWindows{iX}(iQ,2);
                
                iQQ = sum(nQuants(1:iX))+iQ;
                for iR = 1:nRats
                    selectedCells = ratIdx(groupID(:,iA)) == iR;
                    for iE = 1:nEvents
                        tempResponse = nanmean(fullOutputData{iM,iX}{iE,iA}(selectedCells,:),1);
                        meanResponse(iE,iA,iR,iQQ) = nanmean(tempResponse(selectedBins));
                    end
                end
                tempA = reshape(repmat(eventData{iM},1,nRats),[],1);
                tempB = reshape(meanResponse(:,iA,:,iQQ),[],1);
                valid = ~isnan(tempA+tempB);
                [~,~,RL,RU] = corrcoef(tempA(valid),tempB(valid));
                netError(iA,iQQ,1) = RL(2); netError(iA,iQQ,2) = RU(2);
                
                [netCorr(iA,iQQ),netSig(iA,iQQ)] = nancorr(tempA(valid),tempB(valid));
            end
        end
        
%             ShadedErrorbar(1:totalQuants, netCorr(iA,:), [],
%             'L',-netError(iA,:,1)+netCorr(iA,:),'U',netError(iA,:,2)-netCorr(iA,:),'color',GWDColor.(anatGroup{iA}),'lineWidth',1,'marker','none');
        plot(1:totalQuants,netCorr(iA,:),'--o','color',GWDColor.(anatGroup{iA}));
    end
    
    ylim(yRange(iM,:))
    xlim([.5 totalQuants+.5])
    xticks(1:totalQuants);xticklabels(quantNames)
    plotHorizLine(0,{'--k'})
    ylabel('Correlation')
    title(plotTitle{iM})
    
    figFN = [figOutputDir,currSubDir,sprintf('FRto%s_Corr',plotTitle{iM})];
    saveas(fh_data,figFN,'svg'); close(fh_data)
end


