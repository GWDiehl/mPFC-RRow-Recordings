function mPFC_RevisionAnalysis_OfferZone_eLife(varargin)

figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\OfferZone\';

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


% MI for Behavioral Variables

behavioralVar = {'Choice' 'SiteNum' 'OfferDelay' 'OfferValue'};
analysisFN = {'Choice_Accept' 'SiteNumber_ByChoice' 'OfferByChoice_Delay' 'OfferByChoice_Value'};

dataRootName = 'Revision\eLife\MutualInfo\OfferZone\';
dataType = 'MI';

saveSuffix = '_MI';

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
decisionWindow = 4:23;
selectedTime = [8:33];
windowBuffer = [-.5 .5];

selectedEvent = {'OZEntry'};
selectedGroup = {'All'};
plotTitle = {'All'};
nGrps = length(selectedGroup);
fullEventList = repmat(selectedEvent,1,nGrps);

for iV = 1:length(behavioralVar)
    fullDataFN = [dataRootName,analysisFN{iV}];
    
    if strcmp(behavioralVar{iV},'Choice')
        decisionType = {''};
        yRange = [-.015 .24];
        sigYPos = linspace(.15,.22,length(anatGroup))';
        mainDataTitle = {'Choice'};
    else
        if strcmp(behavioralVar{iV},'SiteNum')
            decisionType = {'' '_Accept' '_Skip'};
        else
        decisionType = {'_Accept' '_Skip'};
        end
        yRange = [-.015 .065];   
        sigYPos = linspace(.04,.06,length(anatGroup))';
        mainDataTitle = cellfun(@(x) x(2:end),decisionType,'UniformOutput',0);
    end
    
    nDecision = length(decisionType);
    
    for iD = 1:nDecision
        % Connect to the appropriate decision condition
        currEvent = cellfun(@(x) [x,decisionType{iD}],fullEventList,'UniformOutput',0);
        
        % Full slate of groupings
        [fh_data, outputResults, dataCenters, Params] = PlotFullSingleCellResults(fullDataFN,dataType,currEvent,selectedGroup,plotTitle,anatGroup,...
            'yRange',yRange,'selectedTime',selectedTime,'windowBuffer',windowBuffer);              
        title(mainDataTitle{iD})
        
        % Stats w/ significance for the full data
        [rmANOVA,pVals,sigTime,ANOVA,ANOVA_stats,multiCompare] = ComputeDecisionStats_mPFC(outputResults,anatGroup,decisionWindow);
        
        % Plot the Sig Fig
        [fh_sig] = PlotSingleCellSignificance(sigTime,Params,dataCenters,anatGroup,'xRange',[dataCenters(1)+windowBuffer(1) dataCenters(end)+windowBuffer(2)]);
        
        % Overlay on the original Avg data
        [fh_data] = PlotSingleCellSignificance_OnAvg(fh_data,sigTime,dataCenters,anatGroup,sigYPos);
        
        
        title(mainDataTitle{iD})
        
        set(fh_data,'name',sprintf('MI For %s%s',behavioralVar{iV},decisionType{iD})); 
        set(fh_sig,'name',sprintf('MI For %s%s',behavioralVar{iV},decisionType{iD})); 
        
        figFN = [figOutputDir,currSubDir,sprintf('%s%s_MainData',behavioralVar{iV},decisionType{iD}),saveSuffix];      
        saveas(fh_data,figFN,'svg'); close(fh_data)
        figFN = [figOutputDir,currSubDir,sprintf('%s%s_SigVal',behavioralVar{iV},decisionType{iD}),saveSuffix];     
        saveas(fh_sig,figFN,'svg'); close(fh_sig)
    end    
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Avg Firing rate of cells in each group   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


currSubDir = 'FiringRate\';
SelectDefaultColorMap(parula)

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

anatGroup = {'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
[subA,subB] = computeSubplotSize(length(anatGroup));
climsOverride = {};

selectedTime = [8:33];
decisionBins = [20,33];

selectedEvent = {'OZEntry_Accept' 'OZEntry_Skip'};
selectedGroup = {'All' 'All'};
plotTitle = {'Accept' 'Skip'};
nGrps = length(selectedGroup);

fullDataFN =  'Revision\eLife\PETH\OfferZone\OfferZone_Warped';
dataType = 'General';

[fh_data, outputResults, dataCenters, Params] = PlotFullSingleCellResults(fullDataFN,dataType,selectedEvent,selectedGroup,plotTitle,anatGroup,...
    'selectedTime',selectedTime);
close(fh_data);

transitionBins = cumsum(Params.TimeWindow(1:end-1));
netDiff = cellfun(@(x,y) y-x,outputResults(1,:),outputResults(2,:),'UniformOutput',0);

binsToAvg = dataCenters >= decisionBins(1) & dataCenters <= decisionBins(2);


for iX = 1:3
    fh_data = figure;
    for iA = 1:length(anatGroup)
    subplot(subA,subB,iA)
        [avgAct,order] = sort(nanmean(netDiff{iA}(:,binsToAvg),2));
        
        switch iX
            case 1
                tempData = outputResults{1,iA};
                peak = prctile(tempData(:),[5 95]);
                clims = peak;
                saveSuffix = '_Accept';
            case 2
                tempData = outputResults{2,iA};
                peak = prctile(tempData(:),[5 95]);
                clims = peak;
                saveSuffix = '_Skip';
            case 3                
                tempData = netDiff{iA};
                peak = prctile(abs(avgAct),95);
                clims = [-peak, peak];
                saveSuffix = '_NetSkip';
        end
        s = imagesc(dataCenters,1:length(order),tempData(order,:));
        set(s,'alphadata',~isnan(tempData(order,:)))
        plotVertLine(transitionBins,{'r','linewidth',1.5})
        plotHorizLine(find(avgAct>0,1)-.5,{'--r','linewidth',1.5})
        axis xy
        
        [h,p, chi2stat,df] = prop_test([sum(avgAct>0) round(length(avgAct)/2)] ,repmat(length(avgAct),1,2),0);
        
        if isempty(climsOverride)
            caxis(clims)
        else
            caxis(climsOverride{iA})
        end
        title(anatGroup{iA})
        colorbar
    end
    
    figFN = [figOutputDir,currSubDir,'AvgSortedFR',saveSuffix];
    fh_data.Position = [100 100 1150 620];
    saveas(fh_data,figFN,'svg'); close(fh_data)
end


