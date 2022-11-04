function [behavOutput, FRData, RawFR, MIData, MI_Shuff] = CollectPerSessFRProfiles(inputData,behavData,behavIdx,behavBins,varargin)

nShuff = 0;
spikeBinParams = [5 95 5];
behavBinParams = [5 95 5];

process_varargin(varargin);

%%


nEvents = length(inputData);
nBehav = length(behavData);

[nSess,nLaps,nSites] = size(behavData{1});

FRData = cell(nSess,nBehav);
RawFR = cell(nSess,1);
MIData = cell(nSess,nBehav);
MI_Shuff = cell(nSess,nBehav);
behavOutput = cell(nSess,1);


nCells = cellfun(@(x) size(x,1),inputData{1}(:));

for iS = 1:nSess
    MIData(iS,:) = {nan(nCells(iS),1)};
    MI_Shuff(iS,:) = {nan(nCells(iS),1)};
    
    % Grab the avg FR and behavioral variables for each event sample
    FRProfile = nan(nCells(iS),0);
    behavOutput{iS} = nan(nBehav,0);
    for iE = 1:nEvents
        nSamples = size(inputData{iE}{iS},2);
        FRProfile = cat(2,FRProfile,inputData{iE}{iS});
        tempBehav = nan(nBehav,nSamples);
        for iB = 1:nBehav
            % If all of the behavior is nan just use the event number
            if all(isnan(behavData{iB}(:))) 
                tempBehav(iB,:) = iE;
            else
                temp = reshape(behavData{iB}(iS,:),nLaps,nSites);
                tempBehav(iB,:) = temp(behavIdx{iE}{iS});
            end
        end
        behavOutput{iS} = cat(2,behavOutput{iS},tempBehav);
    end
    
    FRData(iS,:) = arrayfun(@(a) cell2mat(arrayfun(@(x) nanmean(FRProfile(:,behavOutput{iS}(a,:)==x),2),1:behavBins(a),'UniformOutput',0)),1:nBehav,'UniformOutput',0);
    RawFR{iS} = FRProfile;
    
    binnedFR = cell2mat(arrayfun(@(x) discretize(FRProfile(x,:),...
        linspace(prctile(FRProfile(x,:),spikeBinParams(1)),prctile(FRProfile(x,:),spikeBinParams(2)),spikeBinParams(3)+1)),...
        [1:nCells(iS)]','UniformOutput',0));
    
    
    for iB = 1:nBehav
        tempBehav = behavOutput{iS}(iB,:);
        tempBehav = discretize(tempBehav,linspace(prctile(tempBehav,behavBinParams(1)),prctile(tempBehav,behavBinParams(2)),behavBinParams(3)+1));
        for iC = 1:nCells(iS)
            tempData = cat(2,tempBehav',binnedFR(iC,:)');
            MIData{iS,iB}(iC) = calculateInfoMetric(tempData,'PairedMI');
            
            tempOut = nan(nShuff,1);
            for iX = 1:nShuff
                tempData = cat(2,randsample(tempBehav,length(tempBehav))',binnedFR(iC,:)');
                tempOut(iX) = calculateInfoMetric(tempData,'PairedMI');
            end
            MI_Shuff{iS,iB}(iC) = nanmean(tempOut);
        end
    end
    fprintf('Done With %d of %d \n',iS,nSess)
end

MIData = reshape(cell2mat(MIData),sum(nCells),nBehav);
MI_Shuff = reshape(cell2mat(MI_Shuff),sum(nCells),nBehav);

FRData = arrayfun(@(x) cell2mat(FRData(:,x)),1:nBehav,'UniformOutput',0);
