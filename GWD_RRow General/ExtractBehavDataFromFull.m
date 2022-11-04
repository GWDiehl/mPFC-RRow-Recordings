function [dataOfInterest,defaultBinEdges] = ExtractBehavDataFromFull(lapData,sessData,varOfInterest,varargin)

% Extract out a particular variable of interest from the full collected
% lapData/sessData formating

% GWD March 2021

% Dummy Varaible in case you have to pass in specalized info
miscVariable = [];

decisionOptions = {'Skip' 'Earn' 'Quit' 'N/A'};
squeezeD1 = false;

ThresholdZone = 'Global';
ThresholdCalc = 'Heaviside';

process_varargin(varargin);

% If we only passed in a single session worth of data, expand things and
% then collapse them down again at the end
nDim =  ndims(lapData.Decision);
if nDim == 2
    fieldEntries = fieldnames(lapData);
    squeezeD1 = true;
    for iF = 1:length(fieldEntries)
        lapData.(fieldEntries{iF}) = permute(lapData.(fieldEntries{iF}),[3 1 2]);
    end
    if ~isempty(sessData)
        fieldEntries = fieldnames(sessData.Threshold);
        for iF = 1:length(fieldEntries)
            fieldEntries2 = fieldnames(sessData.Threshold.(fieldEntries{iF}));
            for iE = 1:length(fieldEntries2)
                sessData.Threshold.(fieldEntries{iF}).(fieldEntries2{iE}) = permute(sessData.Threshold.(fieldEntries{iF}).(fieldEntries2{iE}),[3 1 2]);
            end
        end
    end
end

WaitZones = find(ismember(lapData.ZoneID(1,1,:),'WaitZone'));
OfferZones = find(ismember(lapData.ZoneID(1,1,:),'OfferZone'));

[LapStart, LapEnd, LapSite, LapAdvance, LapPreStart,...
    OZEnter,OZExit, WZEnter,WZExit, LZEnter,LZExit, AZEnter,AZExit]...
    = identifyRRowLapTimes(lapData,OfferZones,WaitZones,'squeezeD1',0);

[nSess,nLaps,nSites] = size(LapStart);

% If we are wanting past or future lap data
if contains(varOfInterest,'Net')
    temp = varOfInterest;
    varOfInterest = temp(1:(strfind(temp,'Net')-1));
    offset = str2double(temp(strfind(temp,'Net')+3:end));
else
    offset = 0;
end

switch varOfInterest
    case 'Accept'
        dataOfInterest = lapData.AcceptOffer(:,:,OfferZones);
        defaultBinEdges = -.5:1.5;
    case 'Quit'
        dataOfInterest = lapData.QuitOffer(:,:,WaitZones);
        defaultBinEdges = -.5:1.5;
    case 'Quit_Global'
        dataOfInterest = double(lapData.QuitOffer(:,:,WaitZones)==1);
        defaultBinEdges = -.5:1.5;
        
    case 'Skip'
        dataOfInterest = lapData.SkipOffer(:,:,OfferZones);
        defaultBinEdges = -.5:1.5;
    case 'Earn'
        dataOfInterest = lapData.EarnOffer(:,:,WaitZones);
        defaultBinEdges = -.5:1.5;
    case {'Rewarded' 'EarnLap'}
        dataOfInterest = double(lapData.EarnOffer(:,:,WaitZones)==1);
        defaultBinEdges = -.5:1.5;
    case {'UnRewarded' 'NoEarnLap'}
        dataOfInterest = double(lapData.EarnOffer(:,:,WaitZones)~=1);
        defaultBinEdges = -.5:1.5;
        
    case 'EconomicChoice'
        offerValue = ExtractBehavDataFromFull(lapData,[],'Value');
        accept = ExtractBehavDataFromFull(lapData,[],'Accept');
        dataOfInterest = double((offerValue>0 & accept==1) | (offerValue<0 & accept==0));
        defaultBinEdges = -.5:1.5;
    case 'EconomicViolation'
        offerValue = ExtractBehavDataFromFull(lapData,[],'Value');
        accept = ExtractBehavDataFromFull(lapData,[],'Accept');
        dataOfInterest = double((offerValue>0 & accept==0) | (offerValue<0 & accept==1));
        defaultBinEdges = -.5:1.5;
        
        % Accept/Skip/Quit decisions that were adventagous based on the
        % value of the original offer
    case 'GoodAccept'
        goodChoice = ExtractBehavDataFromFull(lapData,[],'EconomicChoice');
        accepts = ExtractBehavDataFromFull(lapData,[],'Accept');
        dataOfInterest = double(goodChoice==1 & accepts==1);
        defaultBinEdges = -.5:1.5;
    case 'GoodSkip'
        goodChoice = ExtractBehavDataFromFull(lapData,[],'EconomicChoice');
        skips = ExtractBehavDataFromFull(lapData,[],'Skip');
        dataOfInterest = double(goodChoice==1 & skips==1);
        defaultBinEdges = -.5:1.5;
    case 'GoodQuit'
        goodChoice = ExtractBehavDataFromFull(lapData,[],'EconomicChoice');
        quits = ExtractBehavDataFromFull(lapData,[],'Quit');
        dataOfInterest = double(goodChoice==1 & quits==1);
        dataOfInterest(isnan(quits)) = nan;
        defaultBinEdges = -.5:1.5;
        
        % Early vs Late quits refer to was the quit made when the remaing
        % delay at the time of quit was above (early) vs below (late)
        % threshold
    case 'EarlyQuit'
        timeToQuit = ExtractBehavDataFromFull(lapData,[],'QuitTime');
        origDelay = ExtractBehavDataFromFull(lapData,[],'Delay');
        delayRemaining = origDelay - timeToQuit;
        threshold = ExtractBehavDataFromFull(lapData,sessData,'Threshold');
        quits = ExtractBehavDataFromFull(lapData,[],'Quit');
        
        dataOfInterest = double(delayRemaining > threshold);
        dataOfInterest(isnan(quits)) = nan;
        defaultBinEdges = -.5:1.5;
    case 'LateQuit'
        timeToQuit = ExtractBehavDataFromFull(lapData,[],'QuitTime');
        origDelay = ExtractBehavDataFromFull(lapData,[],'Delay');
        delayRemaining = origDelay - timeToQuit;
        threshold = ExtractBehavDataFromFull(lapData,sessData,'Threshold');
        quits = ExtractBehavDataFromFull(lapData,[],'Quit');
        
        dataOfInterest = double(delayRemaining < threshold);
        dataOfInterest(isnan(quits)) = nan;
        defaultBinEdges = -.5:1.5;
        
        % Accept/Skip/Quit decisions that were disadventagous based on the
        % value of the original offer
    case 'BadAccept'
        badChoice = ExtractBehavDataFromFull(lapData,[],'EconomicViolation');
        accepts = ExtractBehavDataFromFull(lapData,[],'Accept');
        dataOfInterest = double(badChoice==1 & accepts==1);
        defaultBinEdges = -.5:1.5;
    case 'BadAccept_Earned'
        badChoice = ExtractBehavDataFromFull(lapData,[],'BadAccept');
        earn = ExtractBehavDataFromFull(lapData,[],'Earn');
        dataOfInterest = double(badChoice==1 & earn==1);
        defaultBinEdges = -.5:1.5;
    case 'BadAccept_Unearned'
        badChoice = ExtractBehavDataFromFull(lapData,[],'BadAccept');
        earn = ExtractBehavDataFromFull(lapData,[],'Earn');
        dataOfInterest = double(badChoice==1 & earn==0);
        defaultBinEdges = -.5:1.5;
        
    case 'BadSkip'
        badChoice = ExtractBehavDataFromFull(lapData,[],'EconomicViolation');
        skips = ExtractBehavDataFromFull(lapData,[],'Skip');
        dataOfInterest = double(badChoice==1 & skips==1);
        defaultBinEdges = -.5:1.5;
    case 'BadQuit'
        badChoice = ExtractBehavDataFromFull(lapData,[],'EconomicViolation');
        quits = ExtractBehavDataFromFull(lapData,[],'Quit');
        dataOfInterest = double(badChoice==1 & quits==1);
        dataOfInterest(isnan(quits)) = nan;
        defaultBinEdges = -.5:1.5;
        
    case 'Delay'
        dataOfInterest = lapData.ZoneDelay(:,:,OfferZones);
        defaultBinEdges = 0.5:2:30.5;
    case 'Value'
        dataOfInterest = lapData.Value_H(:,:,OfferZones);
        defaultBinEdges = linspace(-20,20,21);
    case 'Value_Sigmoid'
        dataOfInterest = lapData.Value_H(:,:,OfferZones);
        defaultBinEdges = linspace(-20,20,21);
    case 'Threshold'
        dataOfInterest = permute(repmat(sessData.Threshold.Global.Heaviside(:,WaitZones),1,1,nLaps),[1 3 2]);
        defaultBinEdges = linspace(0,31,32);
    case 'Threshold_Sigmoid'
        dataOfInterest = permute(repmat(sessData.Threshold.Global.Sigmoid(:,WaitZones),1,1,nLaps),[1 3 2]);
        defaultBinEdges = linspace(0,31,32);
    case 'Threshold_SigmoidSlope'
        dataOfInterest = permute(repmat(sessData.Threshold.Global.Sigmoid_Slope(:,WaitZones),1,1,nLaps),[1 3 2]);
        defaultBinEdges = linspace(0,31,32);
        
    case 'GoodDeal'
        offerValue = ExtractBehavDataFromFull(lapData,[],'Value');
        dataOfInterest = double(offerValue>0);
    case 'BadDeal'
        offerValue = ExtractBehavDataFromFull(lapData,[],'Value');
        dataOfInterest = double(offerValue<0);
    case 'OfferEase'
        offerValue = ExtractBehavDataFromFull(lapData,[],'Value');
        dataOfInterest = abs(offerValue);
    case 'OfferDifficulty'
        offerValue = ExtractBehavDataFromFull(lapData,[],'Value');
        dataOfInterest = -abs(offerValue);
        
    case 'ReactionTime'
        dataOfInterest = OZExit - OZEnter;
        defaultBinEdges = linspace(0,4,21);
    case 'SkipRT'
        dataOfInterest = selectData(OZExit - OZEnter,lapData.SkipOffer(:,:,OfferZones)==1);
        defaultBinEdges = linspace(0,4,21);
    case 'AcceptRT'
        dataOfInterest = selectData(OZExit - OZEnter,lapData.AcceptOffer(:,:,OfferZones)==1);
        defaultBinEdges = linspace(0,4,21);
    case 'WZTime'
        dataOfInterest = WZExit - WZEnter;
        defaultBinEdges = linspace(0,30,21);
    case {'LingerTime' 'ConsuptionTime'}
        dataOfInterest = LZExit - LZEnter;
        defaultBinEdges = linspace(3,20,21);
    case {'PrctLingerTime' 'PercentLingerTime'} % Percent of the range of linger times for the given session
        lingerTime = ExtractBehavDataFromFull(lapData,[],'LingerTime');
        dataOfInterest = nan(size(lingerTime));
        for iS = 1:nSess
            temp = lingerTime(iS,:);
            pct = prctile(temp,[1 99]);
            dataOfInterest(iS,:) = (temp - pct(1))/diff(pct)*100;
        end
        defaultBinEdges = linspace(0,100,21);
    case 'QuitTime'
        dataOfInterest = selectData(WZExit - WZEnter,lapData.QuitOffer(:,:,WaitZones)==1);
        defaultBinEdges = linspace(0,20,21);
    case 'AdvanceTime'
        dataOfInterest = AZExit - AZEnter;
        defaultBinEdges = linspace(0,4,21);
    case 'TotalResponseTime'
        dataOfInterest = lapData.TotalSiteTime(:,:,WaitZones);
        defaultBinEdges = 0:1.5:50;
    case 'TotalRestaurantTime'
        dataOfInterest = LapEnd - LapStart;
        defaultBinEdges = 0:1.5:50;
        
    case 'OZRunSpeed'
        dataOfInterest = lapData.RunSpeed(:,:,OfferZones);
        defaultBinEdges = linspace(5,50,21);
    case 'OZPauseTime'
        dataOfInterest = lapData.PauseTime(:,:,OfferZones);
        defaultBinEdges = linspace(0,20,21);
    case 'FullWZRunSpeed'
        dataOfInterest = lapData.RunSpeed(:,:,WaitZones);
        defaultBinEdges = linspace(5,50,21);
    case 'FullWZPauseTime'
        dataOfInterest = lapData.PauseTime(:,:,WaitZones);
        defaultBinEdges = linspace(0,20,21);
        
    case 'WZRunSpeed'
        pathFN = miscVariable;
        dataOfInterest = CollectRRowZoneSpeed(pathFN,WZEnter,WZExit);
        defaultBinEdges = linspace(5,50,21);
    case 'WZPauseTime'
        pathFN = miscVariable;
        [~,dataOfInterest] = CollectRRowZoneSpeed(pathFN,WZEnter,WZExit);
        defaultBinEdges = linspace(0,20,21);
    case 'LZRunSpeed'
        pathFN = miscVariable;
        dataOfInterest = CollectRRowZoneSpeed(pathFN,LZEnter,LZExit);
        defaultBinEdges = linspace(5,50,21);
    case 'LZPauseTime'
        pathFN = miscVariable;
        [~,dataOfInterest] = CollectRRowZoneSpeed(pathFN,LZEnter,LZExit);
        defaultBinEdges = linspace(0,20,21);
    case 'AZRunSpeed'
        pathFN = miscVariable;
        dataOfInterest = CollectRRowZoneSpeed(pathFN,AZEnter,AZExit);
        defaultBinEdges = linspace(5,50,21);
    case 'AZPauseTime'
        pathFN = miscVariable;
        [~,dataOfInterest] = CollectRRowZoneSpeed(pathFN,AZEnter,AZExit);
        defaultBinEdges = linspace(0,20,21);
        
    case 'TaskRunSpeed'
        pathFN = miscVariable;
        dataOfInterest = CollectRRowZoneSpeed(pathFN,LapStart,LapEnd);
        defaultBinEdges = linspace(5,50,21);
    case 'TaskPauseTime'
        pathFN = miscVariable;
        [~,dataOfInterest] = CollectRRowZoneSpeed(pathFN,LapStart,LapEnd);
        defaultBinEdges = linspace(0,20,21);
        
    case {'Decision', 'Choice'}
        choiceStr = lapData.Decision(:,:,:);
        dataOfInterest = cellfun(@(x) find(ismember(decisionOptions,x)),choiceStr);
        dataOfInterest(dataOfInterest == find(ismember(decisionOptions,'N/A'))) = nan;
        defaultBinEdges = 0.5:3.5;
        
        
        % Look at the Wait Zones as those give the actual site number
    case {'Site' 'SiteNumber'}
        dataOfInterest = lapData.ZoneNum(:,:,WaitZones);
        defaultBinEdges = 0.5:4.5;
    case {'Site1' 'Site2' 'Site3' 'Site4'}
        numOfInterest = str2num(varOfInterest(5:end));
        dataOfInterest = lapData.ZoneNum(:,:,WaitZones)==numOfInterest;
        defaultBinEdges = -.5:1.5;
        
        % Ranking of restaurant by threshold
    case {'Rank' 'SiteRank'}
        Threshold = sessData.Threshold.(ThresholdZone).(ThresholdCalc);
        if nDim == 3
            Threshold = permute(repmat(Threshold,1,1,nLaps),[1 3 2]);
        end
        Threshold = Threshold(:,1,:);
        TotalEarns = nansum(lapData.EarnOffer(:,:,WaitZones),2);
        AvgWait = nanmean(selectData(lapData.ZoneDelay(:,:,WaitZones),lapData.EarnOffer(:,:,WaitZones)==1),2);
        
        % Priority: Threshold, Total Rewards, Avg Time waited
        % Ranking data as Sess x Sites x RankVariable
        RankingData = permute(cat(2,Threshold,TotalEarns,AvgWait),[1 3 2]);
        
        [FeederRank,RankOrder] = ExtractRRowRankOrdering(RankingData);
        dataOfInterest = permute(repmat(FeederRank,1,1,nLaps),[1 3 2]);
        defaultBinEdges = 0.5:4.5;
    case {'Rank1' 'Rank2' 'Rank3' 'Rank4'}
        FeederRank = ExtractBehavDataFromFull(lapData,sessData,'Rank');
        
        numOfInterest = str2num(varOfInterest(5:end));
        dataOfInterest = double(permute(repmat(FeederRank,1,1,nLaps),[1 3 2]) == numOfInterest);
        defaultBinEdges = -.5:1.5;
        
    case 'SessionTime'
        sessStart = min(lapData.EnteringZoneTime(:,:),[],2);
        dataOfInterest = lapData.EnteringZoneTime(:,:,OfferZones) - sessStart;
        defaultBinEdges = linspace(0,3600,21);
    case 'LapNumber'
        dataOfInterest = lapData.CurrentLap(:,:,OfferZones);
        defaultBinEdges = linspace(0,250,21);
    case 'RelativeLap'
        dataOfInterest = lapData.CurrentLap(:,:,OfferZones)./max(lapData.CurrentLap(:,:),[],2);
        defaultBinEdges = linspace(0,1,21);
    case 'CycleNumber'
        dataOfInterest = lapData.CurrentCycle(:,:,OfferZones);
        defaultBinEdges = linspace(0,125,21);
    case 'RelativeCycle'
        dataOfInterest = lapData.CurrentCycle(:,:,OfferZones)./max(lapData.CurrentCycle(:,:),[],2);
        defaultBinEdges = linspace(0,1,21);
        
    case {'EarnedRewards' 'TotalEarns'}
        earnLaps = reshape(permute(lapData.EarnOffer(:,:,WaitZones) == 1,[1 3 2]),nSess,[]);
        totalEarns = cumsum(earnLaps,2);
        dataOfInterest = permute(reshape(totalEarns,nSess,nSites,[]),[1 3 2]);
        defaultBinEdges = linspace(0,125,21);
    case 'PercentTotalEarns'
        earns = ExtractBehavDataFromFull(lapData,[],'EarnedRewards');
        maxEarns = max(earns(:,:),[],2);
        dataOfInterest = earns./maxEarns*100;
        defaultBinEdges = linspace(0,100,21);
        
    case 'LapsFromLastEarn'
        earnLaps = reshape(permute(lapData.EarnOffer(:,:,WaitZones) == 1,[1 3 2]),nSess,[]);
        tempData = nan(size(earnLaps));
        for iS = 1:nSess
            validLaps = find(earnLaps(iS,:));
            for iV = 1:length(validLaps)-1
                lapDist = 1:(validLaps(iV+1) - validLaps(iV));
                tempData(iS,validLaps(iV)+1:validLaps(iV+1)) = lapDist;
            end
        end
        dataOfInterest = permute(reshape(tempData,nSess,nSites,[]),[1 3 2]);
        defaultBinEdges = [.5:10.5];
    case 'TimeFromLastEarn'
        earnLaps = reshape(permute(lapData.EarnOffer(:,:,WaitZones) == 1,[1 3 2]),nSess,[]);
        fullLapTimes = reshape(permute(LapStart,[1 3 2]),nSess,[]);
        fullRewardTimes = reshape(permute(LZEnter,[1 3 2]),nSess,[]);
        tempData = nan(size(earnLaps));
        for iS = 1:nSess
            validLaps = find(earnLaps(iS,:));
            for iV = 1:length(validLaps)-1
                lapDist = diff([fullRewardTimes(iS,validLaps(iV)) fullLapTimes(iS,validLaps(iV)+1:validLaps(iV+1))]);
                tempData(iS,validLaps(iV)+1:validLaps(iV+1)) = cumsum(lapDist);
            end
        end
        dataOfInterest = permute(reshape(tempData,nSess,nSites,[]),[1 3 2]);
        defaultBinEdges = linspace(5,100,21);
        
    case 'IdPhi'
        dataOfInterest = lapData.IdPhi(:,:,OfferZones);
    case 'AvgDPhi'
        dataOfInterest = lapData.AvgDPhi(:,:,OfferZones);
        
    case 'PhaseDay'
        dataOfInterest = nan(size(lapData.Decision));
        for iS = 1:size(dataOfInterest,1)
            dataOfInterest(iS,:) = sessData.PhaseDay(iS);
        end
        defaultBinEdges = 0.5:1:14.5;
        
    case 'TrainingDay'
        dataOfInterest = nan(size(lapData.Decision));
        for iS = 1:size(dataOfInterest,1)
            dataOfInterest(iS,:) = sessData.TrainingDay(iS);
        end
        defaultBinEdges = 0.5:1:25.5;
        
    case {'Random' 'RandomLap'}
        dataOfInterest = rand(size(lapData.Decision));
        defaultBinEdges = linspace(0,1,31);
        
    case {'Null' 'All'}
        dataOfInterest = double(~isnan(lapData.ZoneDelay(:,:,OfferZones)));
        defaultBinEdges = [.5 1.5];
        
        
    otherwise
        error('%s Not Implemented',varOfInterest)
end

% Ensure that any laps that were not visited (no delay info) are set to nan
dataOfInterest(isnan(lapData.ZoneDelay(:,:,OfferZones))) = nan;

% If we are looking at data that is forward or backwards in time (Net+:
% Past; Net-: Future) offset the data now
if offset
    dataOfInterest = advanceRRowData(dataOfInterest,offset);
end

if squeezeD1
    dataOfInterest = permute(dataOfInterest,[2 3 1]);
end

