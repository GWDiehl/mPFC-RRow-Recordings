
function [Results, Params] = RRow_RegressionToBehavior(varargin)


% Prelim Directoires & DataDef
[rootDir, DataDefName, BehavDir, TaskSuffix, PathSuffix] = ID_GWD_Dirs;


dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');
subDataDef = subsetFromDataDef(dataDef,OI);

SSNs = subDataDef.SSN;


%%

%%% What kind of analysis do we want to be doing??
% Looking at the avg activity/behavior for a given RRow "Lap" defined by a
% single start/end pairing for the lap?

% Binning data and comparing actiivty/behavior for each bin. Note that
% periods of binned time can be exluded if handled properly.

% Taking a series of arbitrary task defined time chunks and looking at avg
% activity in that chunk (note that behavior should be identical across the
% chunk but there can be multiple chunks for a given behavioral lap

analysisType = 'ByLap'; % 'ByLap' 'ByBin' 'ByChunk'
analysesToRun = {'Correlation' 'SingleRegression' 'MutualInfo' 'StepwiseRegression'};

%%% What window of the data do we want to use for our analysis:

% ByLap: What is the start and end time of each lap (2 element cell)

% ByBin: What is the dt time resolution of our binned spiking data (2
% element cell, first element is the dt, second is passed into
% 'splitRRowTSDs_V2' to restrict the behavioral TSD; use 'All'/'None'/''
% for no restriction)

% ByChunk: What are the start and end times of each chunk we want to
% consider (2xN element cell with start and end points of each chunk)

analysisWindow = {'LapStart' 'LapEnd'};

behavVar = {'SessionTime' 'Value' 'Rank'};
regressionVar = {};
regressOutOrder = [1];

includeData = [];
inInitialModel = [];

dataRegressOrder = 1;


%% General Parameter Values

discretizeData = 0;
MI_Bins = [5 5]; % Max bins for MI calculation for behavior and firing rates
scaleData = 'on';

nShuffles = 0;

individualVar = {'BetaWeight' 'Correlation' 'pVal' 'MutualInfo'};
stepwiseVar = {'BetaWeight' 'pVal' 'InModel'};
stepwiseOrderedVar = {'BetaWeight' 'Correlation' 'pVal'};

process_varargin(varargin);

%%

nSess = length(SSNs);
nBehavVar = length(behavVar);
nRegressVar = length(regressionVar);

if isempty(includeData)
    includeData = true(size(analysisWindow,2),nBehavVar);
end
if isempty(inInitialModel)
    inInitialModel = false(nBehavVar,1);
end


Results = struct;
for iX = 1:length(individualVar)
    Results.Individual.(individualVar{iX}) = cell(nSess,1);
end
for iX = 1:length(stepwiseVar)
    Results.Stepwise.(stepwiseVar{iX}) = cell(nSess,1);
end
for iX = 1:length(stepwiseOrderedVar)
    Results.OrderedStep.(stepwiseOrderedVar{iX}) = cell(nSess,1);
end
if nShuffles
    for iX = 1:length(individualVar)
        Results.Shuffle.Individual.(individualVar{iX}).Example = cell(nSess,1);
        Results.Shuffle.Individual.(individualVar{iX}).Avg = cell(nSess,1);
        Results.Shuffle.Individual.(individualVar{iX}).Std = cell(nSess,1);
    end
    for iX = 1:length(stepwiseVar) 
        Results.Shuffle.Stepwise.(stepwiseVar{iX}).Example = cell(nSess,1);
        Results.Shuffle.Stepwise.(stepwiseVar{iX}).Avg = cell(nSess,1);
        Results.Shuffle.Stepwise.(stepwiseVar{iX}).Std = cell(nSess,1);
    end
    for iX = 1:length(stepwiseOrderedVar)        
        Results.Shuffle.OrderedStep.(stepwiseOrderedVar{iX}).Example = cell(nSess,1);
        Results.Shuffle.OrderedStep.(stepwiseOrderedVar{iX}).Avg = cell(nSess,1);
        Results.Shuffle.OrderedStep.(stepwiseOrderedVar{iX}).Std = cell(nSess,1);
    end
end


% Parameters
Params = struct;
Params.SSN = SSNs';
Params.AnalysisType = analysisType;
Params.AnalysisWindow = analysisWindow;
Params.BehavVar = behavVar;
Params.RegressionVar = regressionVar;
if strcmp(analysisType,'ByChunk')
    Params.IncludedData = includeData;
end

Params.DiscretizeData = discretizeData;
Params.MI_Bins = MI_Bins;
Params.StepwiseZScore = scaleData;

CellID = cell(nSess,1);


%% Get any behavioral Task Data that may be needed

[FullLapData, FullSessionData] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});

[nSess,nLaps,nZones] = size(FullLapData.Decision);
waitZones = find(ismember(squeeze(FullLapData.ZoneID(1,1,:)),'WaitZone'));
offerZones = find(ismember(squeeze(FullLapData.ZoneID(1,1,:)),'OfferZone'));

fullDataDir = cellfun(@(x,y) [rootDir,BehavDir,x,'\',y,'\'],subDataDef.RatID,subDataDef.Directory,'uniformoutput',0);
fullPathFN = cellfun(@(x,y) [x,y,PathSuffix],fullDataDir,SSNs,'UniformOutput',0);

% Grab the relevant Behav Times
[~,~,~,~,~,~,~,~,~,~,~,~,~,eventTime]...
    = identifyRRowLapTimes(FullLapData,offerZones,waitZones);

if strcmp(analysisType,'ByLap') % Extract out the behavioral data for each lap
    fullBehavData = cell(nBehavVar,1);
    fullRegressData = cell(nRegressVar,1);
    for iB = 1:nBehavVar
        if contains(behavVar{iB},'Net')
            temp = behavVar{iB};
            tempVar = temp(1:(strfind(temp,'Net')-1));
            
            [fullBehavData{iB},defaultEdges] = ExtractBehavDataFromFull(FullLapData,FullSessionData,tempVar,'miscVariable',fullPathFN);
            fullBehavData{iB} = advanceRRowData(fullBehavData{iB},str2double(temp(strfind(temp,'Net')+3:end)));
        else
            [fullBehavData{iB},defaultEdges] = ExtractBehavDataFromFull(FullLapData,FullSessionData,behavVar{iB},'miscVariable',fullPathFN);
        end
        if discretizeData
            fullBehavData{iB} = discretize(fullBehavData{iB},defaultEdges);
        end
    end
    for iR = 1:nRegressVar
        if contains(regressionVar{iR},'Net')
            temp = regressionVar{iR};
            tempVar = temp(1:(strfind(temp,'Net')-1));
            
            [fullRegressData{iR},defaultEdges] = ExtractBehavDataFromFull(FullLapData,FullSessionData,tempVar,'miscVariable',fullPathFN);
            fullRegressData{iR} = advanceRRowData(fullRegressData{iR},str2double(temp(strfind(temp,'Net')+3:end)));
        else
            [fullRegressData{iR},defaultEdges] = ExtractBehavDataFromFull(FullLapData,FullSessionData,regressionVar{iR},'miscVariable',fullPathFN);
        end
        if discretizeData
            fullRegressData{iR} = discretize(fullRegressData{iR},defaultEdges);
        end
    end
else
    fullBehavData = cell(0);
    fullRegressData = cell(0);
end

%% Start Looping through each individual behavioral session

OutputsIndividual = cell(nSess,length(individualVar));
OutputsStepwise = cell(nSess,length(stepwiseVar));
OutputOrderedStep = cell(nSess,length(stepwiseOrderedVar));

OutputsIndividual_Shuff = cell(nSess,length(individualVar),3);
OutputsStepwise_Shuff = cell(nSess,length(stepwiseVar),3);
OutputOrderedStep_Shuff = cell(nSess,length(stepwiseOrderedVar),3);

parfor iS = 1:nSess
    [vt, ~, S, S_FNs] = GWD_LoadRRowSession(SSNs{iS},{'Path' 'Spikes'},'restrictToTaskTime',1);
    
    sessStart = vt.x.starttime;
    sessEnd = vt.x.endtime;
    
    nCells = length(S_FNs);
    CellID{iS} = S_FNs;
    
    currLapData = subsetStructByIndx(FullLapData,iS,'squeezeD1',1);
    currSessData = subsetStructByIndx(FullSessionData,iS,'squeezeD1',1);
    lapStart = reshape(eventTime.LapStart(iS,:),nLaps,nZones);
    lapEnd = reshape(eventTime.LapEnd(iS,:),nLaps,nZones);
    
    %%
    
    % Collect up all of the behavioral data and the firing rate data for
    % this session
    
    behavData = cell(nBehavVar,1);
    regressData = cell(nRegressVar,1);
    
    switch analysisType
        case 'ByLap' % Extract out the behavioral data for each lap
            % Just grab the behav data for this session
            behavData = cellfun(@(x) reshape(x(iS,:),nLaps,nZones),fullBehavData,'UniformOutput',0);
            regressData = cellfun(@(x) reshape(x(iS,:),nLaps,nZones),fullRegressData,'UniformOutput',0);
            
            % Compute the avg FR in each lap
            firingRate = nan(nCells,nLaps,nZones);
            currWindowStart = reshape(eventTime.(analysisWindow{1})(iS,:),nLaps,nZones);
            currWindowEnd = reshape(eventTime.(analysisWindow{2})(iS,:),nLaps,nZones);
            windowDur = currWindowEnd - currWindowStart;
            validWindows = ~isnan(windowDur);
            for iC = 1:nCells
                firingRate(iC,validWindows) = arrayfun(@(x,y,z) histcn(S{iC}.range,[x,y])/z,...
                    currWindowStart(validWindows),currWindowEnd(validWindows),windowDur(validWindows));
            end
            binTimes = currWindowStart+windowDur/2;
            timeBinsToInclude = true(numel(binTimes),nBehavVar);
            
        case {'ByBin' 'ByChunk'}
            % Behavior is identical for the two as TSDs, but ByBin can also
            % restrict in various ways (below). Presumibly Chunking times restricts
            % already
            for iP = 1:nBehavVar
                [behavData{iP}, TSDBins] = RRow_TSD_Shell_V2(behavVar{iP},lapStart,lapEnd,...
                    'vt',vt,'LapData',currLapData,'SessData',currSessData,'fullTS',vt.x.range);
                
                if discretizeData
                    behavData{iP}.D = discretize(behavData{iP}.D,linspace(TSDBins(1),TSDBins(2),TSDBins(3)+1));
                end
            end
            for iR = 1:nRegressVar
                [regressData{iR}, TSDBins] = RRow_TSD_Shell_V2(regressionVar{iR},lapStart,lapEnd,...
                    'vt',vt,'LapData',currLapData,'SessData',currSessData,'fullTS',vt.x.range);
                
                if discretizeData
                    regressData{iR}.D = discretize(regressData{iR}.D,linspace(TSDBins(1),TSDBins(2),TSDBins(3)+1));
                end
            end
            switch analysisType
                case 'ByBin'
                    for iP = 1:nBehavVar
                        behavData{iP} = splitRRowTSDs_V2(behavData{iP},analysisWindow{2},lapStart,lapEnd,'LapData',currLapData);
                    end
                    for iR = 1:nRegressVar
                        regressData{iR} = splitRRowTSDs_V2(regressData{iR},analysisWindow{2},lapStart,lapEnd,'LapData',currLapData);
                    end
                    
                    QMtx = MakeQfromS(S,analysisWindow{1},'tStart',sessStart,'tEnd',sessEnd);
                    QMtx = recenterQMtx(QMtx);
                    firingRate = QMtx.data';
                    binTimes = QMtx.range;
                    timeBinsToInclude = true(numel(binTimes),nBehavVar);
                    
                case 'ByChunk'
                    [chunkStart,~,idxID] = IdentifySelectedRRowEvents(currLapData,analysisWindow(1,:));
                    chunkEnd = IdentifySelectedRRowEvents(currLapData,analysisWindow(2,:));
                    chunkDur = chunkEnd - chunkStart;
                    
                    validPairs = ~isnan(chunkDur);
                    nChunks = sum(validPairs(:));
                    firingRate = nan(nCells,nChunks);
                    binTimes = chunkStart(validPairs)+chunkDur(validPairs)/2;
                    idxID = idxID(validPairs);
                    
                    % Logical matrix that signifies for each binTime and
                    % behavior if the data should be used subsequently.
                    % This is originally defined in the input 'includeData'
                    % corresponding to chunk windows and behavioral
                    % variables. Here we expand the general chunk windows
                    % into the specific itterations of each.
                    timeBinsToInclude = false(numel(binTimes),nBehavVar);
                    for iX = 1:size(analysisWindow,2)
                        for iY = 1:nBehavVar
                            timeBinsToInclude(idxID==iX,iY) = includeData(iX,iY);
                        end
                    end
                    
                    for iC = 1:nCells
                        firingRate(iC,:) = arrayfun(@(x,y,z) histcn(S{iC}.range,[x,y])/z,...
                            chunkStart(validPairs),chunkEnd(validPairs),chunkDur(validPairs));
                    end
            end
    end
    
    %% Go run the analysis comparing firing rate to behavior
    
    betaWeight = nan(nCells,nBehavVar,nShuffles+1);
    corrVal = nan(nCells,nBehavVar,nShuffles+1);
    pVal = nan(nCells,nBehavVar,nShuffles+1);
    mInfo = nan(nCells,nBehavVar,nShuffles+1);
    
    stepwiseBeta = nan(nCells,nBehavVar,nShuffles+1);
    stepwisePVal = nan(nCells,nBehavVar,nShuffles+1);
    stepwiseInModel = nan(nCells,nBehavVar,nShuffles+1);
    
    orderedStepBeta = nan(nCells,nBehavVar,dataRegressOrder+1,nShuffles+1);
    orderedStepCorr = nan(nCells,nBehavVar,nShuffles+1);
    orderedStepPVal = nan(nCells,nBehavVar,nShuffles+1);
    
    fullBehav = nan(numel(binTimes),nBehavVar);
    for iP = 1:nBehavVar
        switch analysisType
            case 'ByLap'
                currBD = behavData{iP}(:);
            case {'ByBin' 'ByChunk'}
                currBD = behavData{iP}.data(binTimes);
        end
        fullBehav(:,iP) = currBD;
    end
    fullBehav(~timeBinsToInclude) = nan;
    
    for iC = 1:nCells
        currFR = full(firingRate(iC,:)');       
        
        % Regress out variables of interest in listed order
        for iR = 1:nRegressVar
            switch analysisType
                case 'ByLap'
                    currBD = regressData{iR}(:);
                case {'ByBin' 'ByChunk'}
                    currBD = regressData{iR}.data(binTimes);
            end
            
            validTimes = ~isnan(currFR + currBD);
            currFR(~validTimes) = nan;
            temp = polyfit(currBD(validTimes),currFR(validTimes),regressOutOrder(iR));
            predict = polyval(temp,currBD(validTimes));
            currFR(validTimes) = currFR(validTimes) - predict;
        end
        
        
        for iP = 1:nBehavVar
            currBD = fullBehav(:,iP);
            validTimes = ~isnan(currFR + currBD);
            
            if ismember('Correlation',analysesToRun)
                [corrVal(iC,iP,1),pVal(iC,iP,1)] = nancorr(currBD(validTimes),currFR(validTimes));
                for iX = 1:nShuffles
                    [corrVal(iC,iP,iX+1),pVal(iC,iP,iX+1)] = nancorr(randsample(currBD(validTimes),sum(validTimes)),currFR(validTimes));
                end
            end
            if ismember('SingleRegression',analysesToRun)
                temp = polyfit(currBD(validTimes),currFR(validTimes),dataRegressOrder);
                betaWeight(iC,iP,1) = temp(1);
                for iX = 1:nShuffles
                    temp = polyfit(randsample(currBD(validTimes),sum(validTimes)),currFR(validTimes),dataRegressOrder);
                    betaWeight(iC,iP,iX+1) = temp(1);
                end
            end
            if ismember('MutualInfo',analysesToRun)
                inputData = cat(2,currBD(validTimes),currFR(validTimes));
                for iX = 1:2
                    temp = unique(inputData(:,iX));
                    if length(temp) > MI_Bins(iX)
                        inputData(:,iX) = discretize(inputData(:,iX),linspace(prctile(inputData(:,iX),1),prctile(inputData(:,iX),99),MI_Bins(iX)+1));
                    end
                end
                mInfo(iC,iP,1) = calculateInfoMetric(inputData,'PairedMI');
                for iX = 1:nShuffles
                    mInfo(iC,iP,iX+1) = calculateInfoMetric(cat(2,inputData(:,1),randsample(inputData(:,2),sum(validTimes))),'PairedMI');
                end
            end
        end
        
        if ismember('StepwiseRegression',analysesToRun)
            % Run the stepwise regression across all of the included variables
            validTimes = ~isnan(sum(cat(2,fullBehav,currFR),2));
            [b,se,pval,finalmodel] = stepwisefit(fullBehav(validTimes,:),currFR(validTimes),'display','off','scale',scaleData,'inmodel',inInitialModel);
            
            stepwiseBeta(iC,:,1) = b;
            stepwisePVal(iC,:,1) = pval;
            stepwiseInModel(iC,:,1) = finalmodel;
            
            for iX = 1:nShuffles
                tempBehav = fullBehav;
                for iY = 1:nBehavVar
                    validTimes = ~isnan(tempBehav(:,iY));
                    tempBehav(validTimes,iY) = randsample(tempBehav(validTimes,iY),sum(validTimes));
                end
                validTimes = ~isnan(sum(cat(2,tempBehav,currFR),2));
                [b,se,pval,finalmodel] = stepwisefit(tempBehav(validTimes,:),currFR(validTimes),'display','off','scale',scaleData,'inmodel',inInitialModel);
                
                stepwiseBeta(iC,:,iX+1) = b;
                stepwisePVal(iC,:,iX+1) = pval;
                stepwiseInModel(iC,:,iX+1) = finalmodel;
            end
        end
        if ismember('OrderedStepRegression',analysesToRun)
            [b,r,p,residuals] = OrderedStepwiseRegression(fullBehav,currFR,dataRegressOrder);     
           
            orderedStepBeta(iC,:,:,1) = b;
            orderedStepCorr(iC,:,1) = r;
            orderedStepPVal(iC,:,1) = p;
            
            for iX = 1:nShuffles
                tempBehav = fullBehav;
                for iY = 1:nBehavVar
                    validTimes = ~isnan(tempBehav(:,iY));
                    tempBehav(validTimes,iY) = randsample(tempBehav(validTimes,iY),sum(validTimes));
                end
                
                [b,r,p,residuals] = OrderedStepwiseRegression(tempBehav,currFR,dataRegressOrder);
                
                orderedStepBeta(iC,:,:,iX+1) = b;
                orderedStepCorr(iC,:,iX+1) = r;
                orderedStepPVal(iC,:,iX+1) = p;
            end
        end
    end
    
    OutputsIndividual(iS,:) = cat(1,{betaWeight(:,:,1)},{corrVal(:,:,1)},{pVal(:,:,1)},{mInfo(:,:,1)});
    OutputsStepwise(iS,:) = cat(1,{stepwiseBeta(:,:,1)},{stepwisePVal(:,:,1)},{stepwiseInModel(:,:,1)});
    OutputOrderedStep(iS,:) = cat(1,{orderedStepBeta(:,:,:,1)},{orderedStepCorr(:,:,1)},{orderedStepPVal(:,:,1)});
    
    OutputsIndividual_Shuff(iS,:) = cat(1,{betaWeight(:,:,2)},{corrVal(:,:,2)},{pVal(:,:,2)<0.05},{mInfo(:,:,2)},...
        {nanmean(betaWeight(:,:,2:end),3)},{nanmean(corrVal(:,:,2:end),3)},{nanmean(pVal(:,:,2:end)<0.05,3)},{nanmean(mInfo(:,:,2:end),3)},...
        {nanstd(betaWeight(:,:,2:end),[],3)},{nanstd(corrVal(:,:,2:end),[],3)},{nanstd(pVal(:,:,2:end)<0.05,[],3)},{nanstd(mInfo(:,:,2:end),[],3)});
    OutputsStepwise_Shuff(iS,:) = cat(1,{stepwiseBeta(:,:,2)},{stepwisePVal(:,:,2)<0.05},{stepwiseInModel(:,:,2)},...
        {nanmean(stepwiseBeta(:,:,2:end),3)},{nanmean(stepwisePVal(:,:,2:end)<0.05,3)},{nanmean(stepwiseInModel(:,:,2:end),3)},...
        {nanstd(stepwiseBeta(:,:,2:end),[],3)},{nanstd(stepwisePVal(:,:,2:end)<0.05,[],3)},{nanstd(stepwiseInModel(:,:,2:end),[],3)});
    OutputOrderedStep_Shuff(iS,:) = cat(1,{orderedStepBeta(:,:,:,2)},{orderedStepCorr(:,:,2)},{orderedStepPVal(:,:,2)<0.05},...
        {nanmean(orderedStepBeta(:,:,:,2:end),4)},{nanmean(orderedStepCorr(:,:,2:end),3)},{nanmean(orderedStepPVal(:,:,2:end)<0.05,3)},...
        {nanstd(orderedStepBeta(:,:,:,2:end),[],3)},{nanstd(orderedStepCorr(:,:,2:end),[],3)},{nanstd(orderedStepPVal(:,:,2:end)<0.05,[],3)});
    
    fprintf('Done with %s: %d of %d \n',SSNs{iS},iS,nSess)
    
end

Params.CellID = CellID;

% Package up the results
for iX = 1:length(individualVar)
    Results.Individual.(individualVar{iX}) = OutputsIndividual(:,iX);
end
for iX = 1:length(stepwiseVar)
    Results.Stepwise.(stepwiseVar{iX}) = OutputsStepwise(:,iX);
end
for iX = 1:length(stepwiseOrderedVar)
    Results.OrderedStep.(stepwiseOrderedVar{iX}) = OutputOrderedStep(:,iX);
end
if nShuffles
    for iX = 1:length(individualVar)
        Results.Shuffle.Individual.(individualVar{iX}).Example = OutputsIndividual_Shuff(:,iX,1);
        Results.Shuffle.Individual.(individualVar{iX}).Avg = OutputsIndividual_Shuff(:,iX,2);
        Results.Shuffle.Individual.(individualVar{iX}).Std = OutputsIndividual_Shuff(:,iX,3);
    end
    for iX = 1:length(stepwiseVar)
        Results.Shuffle.Stepwise.(stepwiseVar{iX}).Example = OutputsStepwise_Shuff(:,iX,1);
        Results.Shuffle.Stepwise.(stepwiseVar{iX}).Avg = OutputsStepwise_Shuff(:,iX,2);
        Results.Shuffle.Stepwise.(stepwiseVar{iX}).Std = OutputsStepwise_Shuff(:,iX,3);
    end
    for iX = 1:length(stepwiseOrderedVar)
        Results.Shuffle.OrderedStep.(stepwiseOrderedVar{iX}).Example = OutputOrderedStep_Shuff(:,iX,1);
        Results.Shuffle.OrderedStep.(stepwiseOrderedVar{iX}).Avg = OutputOrderedStep_Shuff(:,iX,2);
        Results.Shuffle.OrderedStep.(stepwiseOrderedVar{iX}).Std = OutputOrderedStep_Shuff(:,iX,3);
    end
end
