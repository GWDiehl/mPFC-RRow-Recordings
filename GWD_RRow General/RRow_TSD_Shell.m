function [TSD,TCBins] = RRow_TSD_Shell(TSD_ToMake,LapStart,LapEnd,varargin)

% Shell function that does any necessary first pass processing/prep work
% for building a wide variety of TSD objects for data in the RRow task. In
% general this should be the user call in making TSDs. Builds the
% respecting TSD object and also returns the default binning suggestions
% for the particular variable. Outputs can be passed right to TuningCurve
% or Decoding functions.
%
% Inputs
% TSD_ToMake - Cell Array of string of the TSD that you would like to make.
%       This can be any number of concurrent TSDs or in the case of a
%       single can be input as a simple str.
%       NOTE: Any TSD can be generated with shuffled data entries by adding
%       a "_Shuff" suffix to its str.
% LapStart/LapEnd - A Laps x Sites matrix of the lap start/end times and
%       used for any lap based TSDs
%
% Other Inputs
% ThresholdZone/Calc - Strings that correspond to the fields in LapData
%       that are used for value calcs.
% fullTS - Vector of times of TS/TSD of a desired expansion of the TSD
% bins - Default bins to output for subsetquent TC/Decoding. Some TSDs
%       (categorical data) will override this.
%
% vt - Standard session tracking data required for some TSDs
% LapData/SessData/IdPhi - Single session sturcutres with Lap x Site info
%       relevant for various TSDs (Choice, Times, etc.). THis origianlly
%       comes from 'collectRRowDataToAnalyze' and can be refined down to a
%       sinlge session via 'subsetStructByIndx'
%
% Outputs
% TSD - Cell array of the requested TSD matching the TSD_ToMake input
% TCBins - Cell array with each entry a 3 element vector composing the
%       default recomended [min max nBins] for that TSD

% GWD Feb 2020


%% Currently Supported variables:

% Requires: vt
%%% XSpace, YSpace, Radius, Angle, Velocity, Accelateration, AnlgeRelative,
%%% TimeFromStart, TimeFromEnd, Random, RandomSmooth, SessionTime

% Requires: vt, LapData
%%% RadiusRelative, LapTime, DelayElapsed, PercentDelayElapsed,
%%% DelayRemaining, PercentDelayRemaining

% Requires: LapData
%%% Delay, Skip, Accept, Earn, Quit, NoEarnLap, EarnLap, Choice, MazeChunk,
%%% TotalEarns, PercentTotalEarns

% Requires: LapData, SessData
%%% Value, GoodDeal, BadDeal, SiteRank, NextRank, OfferDifficulty,
%%% OfferEase

% Requires: IdPhi
%%% IdPhi, zIdPhi, AvgDPhi, VTE

% Requires: IdPhi, LapData
%%% VTE_BySite

% Requires: Nothing
%%% SiteNumber, NextSite


%%

% Which threhold to use when calculating value
ThresholdZone = 'Global';
ThresholdCalc = 'Heaviside';

fullTS = [];
bins = 30;

vt = [];

LapData = [];
SessData = [];

process_varargin(varargin);

%%

if ~iscell(TSD_ToMake)
    TSD_ToMake = {TSD_ToMake};
    returnCell = 0;
else
    returnCell = 1;
end

nZones = size(LapStart,2);
WaitZones = 1:nZones;
OfferZones = WaitZones+nZones;

validLaps = ~isnan(LapStart+LapEnd);

if isempty(LapData)    
    SiteNum = nan(size(LapStart));
    LapAdvance = [];
else
    SiteNum = LapData.ZoneNum(:,WaitZones);
    
    AdvanceAccept = LapData.ExitZoneTime(:,WaitZones); % Exits from Accepted Waits (Earn/Quit)
    AdvanceSkip = LapData.ExitZoneTime(:,OfferZones); % Exits from OfferZones
    AdvanceSkip(LapData.AcceptOffer(:,OfferZones)==1) = nan; % Get rid of OfferZone exits when accepting
    tmp = cat(3,AdvanceAccept,AdvanceSkip);
    LapAdvance = squeeze(nanmean(tmp,3)); %Nanmean to take only the real value
end

TCDim = length(TSD_ToMake);
TCBins = cell(1,TCDim);
TSD = cell(1,TCDim);
if length(bins) == 1 && TCDim > 1
    bins = repmat(bins,1,TCDim);
end

restrictToValid = 1;
        
for iD = 1:TCDim
    scratchTSD = 1;
    shuffleBehavData = 0;
    CurrentTC_Var = TSD_ToMake{iD};
    
    StartTimes = LapStart;
    EndTimes = LapEnd;
    SiteNum_In = SiteNum;
    
    if contains(CurrentTC_Var,'_Shuff')
        CurrentTC_Var = strrep(CurrentTC_Var,'_Shuff','');
        shuffleBehavData = 1;
    end
    
    
    switch CurrentTC_Var
        case {'IdPhi' 'zIdPhi' 'VTE' 'VTE_BySite'}
            TSDData = LapData.IdPhi(:,OfferZones);
        case 'AvgDPhi'
            TSDData = LapData.AvgDPhi(:,OfferZones);
        case 'Delay'
            TSDData = LapData.ZoneDelay(:,OfferZones);
        case {'Value' 'GoodDeal' 'BadDeal' 'OfferDifficulty' 'OfferEase'}
            Delay = LapData.ZoneDelay(:,OfferZones);
            TSDData = SessData.Threshold.(ThresholdZone).(ThresholdCalc) - Delay;
            if strcmp(CurrentTC_Var,'GoodDeal') % Positive Value
                TSDData = TSDData >= 0;
            elseif strcmp(CurrentTC_Var,'BadDeal') % Negative Value
                TSDData = TSDData < 0;
            end
        case 'AngleRelative'
            TSDData = LapData.FeederTimes(:,WaitZones);
        case 'RadiusRelative'
%             TSDData = LapData.AcceptOffer(:,OfferZones); Old method
            TSDData = LapData.FeederTimes(:,WaitZones);
        case {'LapTime','LapTime_Earn','LapTime_NoEarn'}
            LapMid = (LapEnd+LapStart)/2; %Middle of the lap
            TSDData = LapMid;
            FeederTimes = LapData.FeederTimes(:,WaitZones);
            RewardLaps = LapData.EarnOffer(:,WaitZones) ==1;
            TSDData(RewardLaps) = FeederTimes(RewardLaps);
            switch CurrentTC_Var
                case 'LapTime_Earn'
                    TSDData(~RewardLaps) = nan;
                case 'LapTime_NoEarn'
                    TSDData(RewardLaps) = nan;
            end
        case {'TimeFromStart_Earn' 'TimeFromEnd_Earn'}
            TSDData = LapData.EarnOffer(:,WaitZones) ==1;
        case {'TimeFromStart_NoEarn' 'TimeFromEnd_NoEarn'}
            TSDData = LapData.EarnOffer(:,WaitZones) ==0;
            
            
        case 'ReactionTime'
            TSDData = LapData.ExitZoneTime(:,OfferZones) - LapData.EnteringZoneTime(:,OfferZones);
            
        case 'ConsumptionTime'
            TSDData = LapAdvance - LapData.FeederTimes(:,WaitZones);
            
        case 'QuitTime'
            TSDData = LapData.ExitZoneTime(:,WaitZones) - LapData.EnteringZoneTime(:,WaitZones);
            TSDData = selectData(TSDData,LapData.QuitOffer(:,WaitZones) == 1);
            
        case {'DelayElapsed' 'DelayRemaining' 'PercentDelayElapsed' 'PercentDelayRemaining'} % Only valid during the wait zone
            StartTimes = LapData.EnteringZoneTime(:,WaitZones);
            EndTimes = min(LapData.ExitZoneTime(:,WaitZones),LapData.FeederTimes(:,WaitZones));
            TSDData = LapData.ZoneDelay(:,WaitZones);
            
        case 'TimeToQuit'
            StartTimes = LapData.EnteringZoneTime(:,WaitZones);
            EndTimes = LapData.ExitZoneTime(:,WaitZones);            
            TSDData = LapData.ExitZoneTime(:,WaitZones) - LapData.EnteringZoneTime(:,WaitZones);
            
            StartTimes = selectData(StartTimes,LapData.QuitOffer(:,WaitZones) == 1);
            EndTimes = selectData(EndTimes,LapData.QuitOffer(:,WaitZones) == 1);
            TSDData = selectData(TSDData,LapData.QuitOffer(:,WaitZones) == 1);
            
        case {'TimeToLeaveFood' 'TimeLingered'}
            StartTimes = LapData.FeederTimes(:,WaitZones);
            EndTimes = LapData.ExitZoneTime(:,WaitZones);            
            TSDData = LapAdvance - LapData.FeederTimes(:,WaitZones);
            
            EndTimes = selectData(EndTimes,LapData.EarnOffer(:,WaitZones) == 1);
            
        case 'Choice'
            choiceOptions = {'Skip' 'Earn' 'Quit'};
            lapChoices = LapData.Decision;
            TSDData = nan(size(lapChoices));
            for iX = 1:length(choiceOptions)
                TSDData(ismember(lapChoices,choiceOptions{iX})) = iX;
            end   
            
        case 'NextSite'
            SiteNum_In = advanceRRowData(SiteNum_In,-1);
            TSDData = nan(size(validLaps));
        case 'PreviousSite'
            SiteNum_In = advanceRRowData(SiteNum_In,1);
            TSDData = nan(size(validLaps));
            
        case {'SiteRank' 'NextRank' 'PreviousRank'}
            Threshold = SessData.Threshold.(ThresholdZone).(ThresholdCalc);
            TotalEarns = nansum(LapData.EarnOffer(:,WaitZones));
            AvgWait = nanmean(selectData(LapData.ZoneDelay(:,WaitZones),LapData.EarnOffer(:,WaitZones)==1));
            
            % Priority: Threshold, Total Rewards, Avg Time waited
            RankingData = cat(1,Threshold,TotalEarns,AvgWait);
            
            [~,FeederRank] = sortrows(RankingData','descend');
            origSite = SiteNum_In;
            SiteNum_In = nan(size(origSite));
            for iS = 1:length(FeederRank)
                SiteNum_In(origSite == FeederRank(iS)) = iS;
            end
            if strcmp(CurrentTC_Var,'NextRank')
                SiteNum_In = advanceRRowData(SiteNum_In,-1);                
            elseif strcmp(CurrentTC_Var,'PreviousRank')
                SiteNum_In = advanceRRowData(SiteNum_In,1);
            end
            TSDData = nan(size(validLaps));
            
            
        case 'MazeChunk'
            chunkOptions = {'OfferZone' 'WaitZone' 'LingerZone' 'AdvanceZone'};
            TSDData = [];
            StartTimes = [];
            restrictToValid = 0;
            for iX = 1:length(chunkOptions)
                switch chunkOptions{iX}
                    case 'OfferZone'
                        tempStart = reshape(LapData.EnteringZoneTime(:,OfferZones),[],1);
                    case 'WaitZone'
                        tempStart = reshape(LapData.EnteringZoneTime(:,WaitZones),[],1);
                    case 'LingerZone'
                        tempStart = reshape(LapData.FeederTimes(:,WaitZones),[],1);
                    case 'AdvanceZone'
                        tempStart = reshape(LapAdvance,[],1);
                end
                tempStart = tempStart(~isnan(tempStart));
                tempData = repmat(iX,size(tempStart));
                
                StartTimes = cat(1,StartTimes,tempStart);
                TSDData = cat(1,TSDData,tempData);
            end
                
        case 'Skip'
            TSDData = LapData.SkipOffer(:,OfferZones);
        case 'Accept'
            TSDData = LapData.AcceptOffer(:,OfferZones);
        case 'PreviousAccept'
            TSDData = advanceRRowData(LapData.AcceptOffer(:,OfferZones),1);
        case 'NextAccept'
            TSDData = advanceRRowData(LapData.AcceptOffer(:,OfferZones),-1);
        case 'Quit'
            TSDData = LapData.QuitOffer(:,WaitZones);
        case 'Earn'
            TSDData = LapData.EarnOffer(:,WaitZones);
            
        case 'NoEarnLap'
            TSDData = LapData.QuitOffer(:,WaitZones)==1 | LapData.SkipOffer(:,OfferZones)==1;
        case 'EarnLap'
            TSDData = LapData.EarnOffer(:,WaitZones) == 1;
        case {'TotalEarns' 'PercentTotalEarns'}
            EarnLaps = LapData.EarnOffer(:,WaitZones) == 1;
            SummedLaps = cumsum(reshape(EarnLaps',[],1));
            TSDData = reshape(SummedLaps,size(EarnLaps,2),size(EarnLaps,1))';
            
        otherwise
            TSDData = nan(size(validLaps));
            scratchTSD = 0; % Data was derived from a TSD
    end
    
    if restrictToValid
        TSDData = TSDData(validLaps);
        StartTimes = StartTimes(validLaps);
        EndTimes = EndTimes(validLaps);
        SiteNum_In = SiteNum_In(validLaps);
    end
    
    if shuffleBehavData && scratchTSD % Before you make the TSD
        TSDData = randsample(TSDData,length(TSDData));
    end
    
    [TSD{iD}, binEdges, nBinsOverRide] = createRRowTSD_GWD(CurrentTC_Var,vt,StartTimes,EndTimes,TSDData,SiteNum_In,'fullTS',fullTS);
    if ~isempty(nBinsOverRide)
        bins(iD) = nBinsOverRide;
    end
    TCBins{iD} = [binEdges bins(iD)];
    
    if shuffleBehavData && ~scratchTSD % Single extracted TSD
        TSD{iD}.D = randsample(TSD{iD}.D,length(TSD{iD}.D));
    end
end
if ~returnCell
    TSD = TSD{:};
    TCBins = TCBins{:};
end
