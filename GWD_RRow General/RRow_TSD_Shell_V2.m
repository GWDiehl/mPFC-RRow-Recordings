function [TSD,TCBins] = RRow_TSD_Shell_V2(TSD_ToMake,LapStart,LapEnd,varargin)

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


% Updated from original function such that per lap data is extracted from
% "ExtractbehavDataFromFull"
% GWD April 2022


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
    
    % If we are wanting past or future lap data
    if contains(CurrentTC_Var,'Net')
        temp = CurrentTC_Var;
        CurrentTC_Var = temp(1:(strfind(temp,'Net')-1));
        offset = str2double(temp(strfind(temp,'Net')+3:end));
    else
        offset = 0;
    end
    
    
    switch CurrentTC_Var
        
        case {}
            EarnLaps = LapData.EarnOffer(:,WaitZones) == 1;
            SummedLaps = cumsum(reshape(EarnLaps',[],1));
            TSDData = reshape(SummedLaps,size(EarnLaps,2),size(EarnLaps,1))';
        
        % These are per lap data that can just be extracted directly
        case {'Skip' 'Accept' 'Quit' 'Earn' 'IdPhi' 'AvgDPhi' 'Delay' 'Value' 'GoodDeal' 'BadDeal' 'ReactionTime' ...
                'ConsumptionTime' 'QuitTime' 'Choice' 'Decision' 'Site' ...
                'SiteNumber' 'PercentTotalEarns' 'TotalEarns' 'LapsFromLastEarn' 'EarnLap' 'NoEarnLap' ...
                'OZRunSpeed' 'EconomicChoice' 'LapsFromLastEarn' 'LingerTime' 'PrctLingerTime'}
            TSDData = ExtractBehavDataFromFull(LapData,SessData,CurrentTC_Var);
            
        case {'Rank' 'SiteRank'}
            TSDData = ExtractBehavDataFromFull(LapData,SessData,CurrentTC_Var,...
                'ThresholdZone',ThresholdZone,'ThresholdCalc',ThresholdCalc);            
            
        case 'RandomLap'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'Random');
            
            % These are per lap data but are further modified before
            % creating the TSD
        case {'zIdPhi' 'VTE' 'VTE_BySite'}
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'IdPhi');
            
        case {'OfferDifficulty' 'OfferEase'}
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'Value');
            
        case {'TimeFromStart_Earn' 'TimeFromEnd_Earn'}
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'Rewarded');            
        case {'TimeFromStart_NoEarn' 'TimeFromEnd_NoEarn'}
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'UnRewarded');
            
            % Make sure that old terminology (Next/Previous) will still work 
        case 'NextSite'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'SiteNet-1');
        case 'PreviousSite'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'SiteNet1'); 
        case 'NextAccept'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'AcceptNet-1');
        case 'PreviousAccept'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'AcceptNet1');
            
        case 'NextRank'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'RankNet-1',...
                'ThresholdZone',ThresholdZone,'ThresholdCalc',ThresholdCalc);
        case 'PreviousRank'
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'RankNet1',...
                'ThresholdZone',ThresholdZone,'ThresholdCalc',ThresholdCalc);
            
            
            % These are things that are distinct for each time bin
        case {'AngleRelative' 'RadiusRelative'}
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
            
        case {'DelayElapsed' 'DelayRemaining' 'PercentDelayElapsed' 'PercentDelayRemaining'} % Only valid during the wait zone
            StartTimes = LapData.EnteringZoneTime(:,WaitZones);
            EndTimes = min(LapData.ExitZoneTime(:,WaitZones),LapData.FeederTimes(:,WaitZones));
            TSDData = ExtractBehavDataFromFull(LapData,SessData,'Delay');
            
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
            
        otherwise
            TSDData = nan(size(validLaps));
            scratchTSD = 0; % Data was derived from a TSD
    end
    
    % If we are looking at per lap data and are future/past
    if offset
        TSDData = advanceRRowData(TSDData,offset);
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
