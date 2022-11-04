function [outputTSD, binEdges, nBins] = createRRowTSD_GWD(Variable,vt,LapStart,LapEnd,LapData,LapSite,varargin)

% Inner working function for construcitng TSDs on the RRow task across a
% wide variety of dimentions. Each call will produce a single TSD along
% with the default set bin edges and number of bins.
%
% TSDs will default to the minimum number of time/data bins to account for
% the variable domain. At a minimum this is a point at the beginning of an
% entry and immediatly before a change to the data field.
%
% Input
% Variable - Str of which TSD you want
% vt - raw tracking data that may be needed for some TSDs
% LapStart - Start times of each lap in the RRow Task, taken as entry into
%       the offer zone at each site. Used for calculating things based on
%       each site.
% LapEnd - End time of each lap in RRow, taken as the entry into the offer
%       zone of the subseqnet site.
% LapSite - RR number of the respective Start/End times
%
% Optional Inputs
% fullTS - A string of TS that can be entred if you would like your output
%       TSD to not be the minimum but rather have entries for a specific
%       set of times. In general this is used to build a high time
%       resolution TSD such as based on the fullTS of the orig tracking
%       data.
% 
% OUTPUTS 
% outputTSD - A generated TSD
% binEdges - The default min/max bin edges for the particualr variable
%       based on empical assesment of a good range of values.
% nBins - An infrequent output only when the constructed TSD is of a
%       categorical nature. Effectivly how many categories exist for the
%       particular varible and what should be used whenever binning based
%       on TSD.

% GWD 2019


%% Currently Supported variables:

%%% XSpace, YSpace, Radius, Angle, Velocity, Accelateration, AnlgeRelative,
%%% TimeFromStart, TimeFromEnd, Random, RandomSmooth, SessionTime

%%% RadiusRelative, LapTime 

%%% Delay, Skip, Accept, Earn, Quit, NoEarnLap, EarnLap, Choice, MazeChunk

%%% Value, GoodDeal, BadDeal, SiteRank, NextRank

%%% IdPhi, zIdPhi, AvgDPhi, VTE, VTE_BySite

%%% SiteNumber, NextSite

%%

fullTS = [];

process_varargin(varargin);

%%

nBins = [];


switch Variable
    % These already exist as TSDs so pull directly
    case 'XSpace'
        outputTSD = vt.x;
        binEdges = [-80 80];
    case 'YSpace'
        outputTSD = vt.y;
        binEdges = [-80 80];
    case 'Radius'
        outputTSD = vt.r;
        binEdges = [30 100];
    case 'Angle'
        outputTSD = vt.p;
        binEdges = [-pi pi];
    case 'Velocity'
        outputTSD = vt.v;
        binEdges = [0 80];
    case 'Acceleration'
        outputTSD = vt.a;
        binEdges = [-100 100];
        
    case 'AngleRelative'
        orig_state = warning;
       warning('off','all')
        % Old method normalizing each lab individually
        % % %          AngleRelative = vt.p;
        % % %         Times = [0; sort(LapStart)];
        % % %         StartAngles = [nan; vt.p.data(sort(LapStart))];
        % % %         EndAngles = [nan; vt.p.data(sort(LapEnd))];
        % % %         AngleRange = circ_dist(EndAngles,StartAngles);
        % % %
        % % %         RelAngle = arrayfun(@(x,y) circ_dist(y,StartAngles(find(Times <= x,1,'last')))...
        % % %             /AngleRange(find(Times <= x,1,'last')),vt.p.range,vt.p.D);
        % % %         AngleRelative.D = RelAngle;
        
        % New method normalizing based on a standard/avged start/end point
        zoneIDs = unique(LapSite);
        AngleRelative = vt.p;
        AngleRelative.D(:) = nan;
        
        startAngles = nan(1,length(zoneIDs));        
        endAngles = nan(1,length(zoneIDs));
        
        for iZ = 1:length(zoneIDs)
            LapsToUse = LapSite == zoneIDs(iZ);
            enterTimes = LapStart(LapsToUse);
            exitTimes = LapEnd(LapsToUse);
            
            temp = vt.p.data(enterTimes);
            startAngles(iZ) = circ_median(temp(~isnan(temp)));
            
            temp = vt.p.data(exitTimes);
            endAngles(iZ) = circ_median(temp(~isnan(temp)));
        end
        angleRange = circ_dist(endAngles,startAngles);
        
        for iL = 1:length(LapStart)
             validIdx = vt.p.range >= LapStart(iL) & vt.p.range <= LapEnd(iL);            
            tempData = vt.p.D(validIdx);            
            scaledData = circ_dist(tempData,startAngles(LapSite(iL)))/angleRange(LapSite(iL));  
            AngleRelative.D(validIdx) = scaledData;
        end
        % Could try to come up with some mechanism to align/normalize the
        % location of the feeders
% % %         rewardAngles = nan(1,length(zoneIDs));
% % %         for iZ = 1:length(zoneIDs)
% % %             LapsToUse = LapSite == zoneIDs(iZ);
% % %             rewardAngles(iZ) = findRewardFromPath({AngleRelative},LapData(LapsToUse));
% % %         end
% % %         
% % %         for iL = 1:length(LapStart)
% % %             tempData = AngleRelative.restrict(LapStart(iL),LapEnd(iL));
% % %             TransitionIdx = find(tempData.D(4:end) > rewardAngles(LapSite(iZ)),1)+3;
% % %             
% % %             preData = tempData.D(1:TransitionIdx)/tempData.D(TransitionIdx)*rewardAngle;
% % %             tempData.D(1:TransitionIdx) = preData;
% % %             
% % %             postData = tempData.D(TransitionIdx+1:end) - tempData.D(TransitionIdx+1);
% % %             postData = postData/postData(end)*(1-rewardAngle)+rewardAngle;
% % %             
% % %             tempData.D(TransitionIdx+1:end) = postData;
% % %             
% % %             AngleRelative.D(ismember(AngleRelative.range, tempData.range)) = tempData.data;
% % %             
% % %         end
        
        outputTSD = AngleRelative;
        binEdges = [-.1 1.1];        
        warning(orig_state)
        
    case 'RadiusRelative'
        zoneIDs = unique(LapSite);
        RadiusRelative = vt.r;
        RadiusRelative.D(:) = nan;
        minRadius = nan(1,length(zoneIDs));
        RadiusRange = nan(1,length(zoneIDs));
        
        for iZ = 1:length(zoneIDs)
            % Old method taking the farthest/closest out in a site
            % % %             % LapData == Accepts; Take accepts at the site of interest
            % % %             LapsToUse = LapSite == zoneIDs(iZ) & LapData == 1;
            % % %             enterTimes = LapStart(LapsToUse);
            % % %             exitTimes = LapEnd(LapsToUse);
            % % %             minRadius(iZ) = nanmedian(arrayfun(@(x,y) prctile(vt.r.restrict(x,y).D,1),enterTimes,exitTimes));
            % % %             maxRadius = nanmedian(arrayfun(@(x,y) prctile(vt.r.restrict(x,y).D,99),enterTimes,exitTimes));
            % % %             RadiusRange(iZ) = maxRadius - minRadius(iZ);
            
            % New method, taking the location at the reward delivery == r1 and at the
            % zone entry == r0
            % LapData == FeederTimes; Take earns at the site of interest
            LapsToUse = LapSite == zoneIDs(iZ) & ~isnan(LapData);
            enterTimes = LapStart(LapsToUse);
            minRadius(iZ) = nanmedian(vt.r.data(enterTimes));
            
            maxRadius = findRewardFromPath({vt.r},LapData(LapsToUse));
            RadiusRange(iZ) = maxRadius - minRadius(iZ);
        end
        
        for iL = 1:length(LapStart)
            validIdx = vt.r.range >= LapStart(iL) & vt.r.range <= LapEnd(iL);            
            tempData = vt.r.D(validIdx);
            scaledData = (tempData - minRadius(LapSite(iL)))/RadiusRange(LapSite(iL));
            RadiusRelative.D(validIdx) = scaledData;
        end
        
        outputTSD = RadiusRelative;
        binEdges = [-.25 1.25];
        
    case {'LapTime','LapTime_Earn','LapTime_NoEarn'} % Time +- from the center of the lap, either reward or midTime (skips)
        CenterPoint = LapData;
        nLaps = length(LapStart);
        LapTime = vt.x;
        LapTime.D(:) = nan;
        timeData = LapTime.range;
        for iL = 1:nLaps
            validIdx = timeData >= LapStart(iL) & timeData <= LapEnd(iL);
            tempTime = timeData(validIdx);
            centeredLap = tempTime - CenterPoint(iL);
            LapTime.D(validIdx) = centeredLap;
        end
        outputTSD = LapTime;
        binEdges = [-25 20];        
        
    case {'TimeFromStart' 'TimeFromEnd' ...
            'TimeFromStart_Earn' 'TimeFromStart_NoEarn'...
            'TimeFromEnd_Earn' 'TimeFromEnd_NoEarn'}
        nLaps = length(LapStart);
        LapTime = vt.x;
        LapTime.D(:) = nan;
        timeData = LapTime.range;
        switch Variable
            case 'TimeFromStart'
                comparePt = LapStart;
                binEdges = [0 7.5];
            case 'TimeFromEnd'
                comparePt = LapEnd;
                binEdges = [-7.5 0];                
        end
        if ismember(Variable,{'TimeFromStart_Earn' 'TimeFromStart_NoEarn'...
            'TimeFromEnd_Earn' 'TimeFromEnd_NoEarn'})
            comparePt(~LapData) = nan;
        end
        for iL = 1:nLaps
            validIdx = timeData >= LapStart(iL) & timeData <= LapEnd(iL);
            tempTime = timeData(validIdx);
            currTime = tempTime - comparePt(iL);
            LapTime.D(validIdx) = currTime;
        end
        outputTSD = LapTime;
        
    case {'DelayElapsed' 'PercentDelayElapsed' 'TimeLingered'}
        nLaps = length(LapStart);
        LapTime = vt.x;
        LapTime.D(:) = nan;
        timeData = LapTime.range;
        for iL = 1:nLaps
            if isnan(LapStart(iL)) % Delay fields will have nan for skip laps
                continue
            end
            validIdx = timeData >= LapStart(iL) & timeData <= LapEnd(iL);
            tempTime = timeData(validIdx);
            currTime = tempTime - LapStart(iL);
            if strcmp(Variable,'PercentDelayElapsed')
                currTime = currTime/LapData(iL)*100;
            end
            LapTime.D(validIdx) = currTime;
        end
        outputTSD = LapTime; 
        switch Variable
            case 'PercentDelayElapsed'
                binEdges = [0 100];
            case 'DelayElapsed'
                binEdges = [0 30];
            case 'TimeLingered'
                binEdges = [0 34];
        end
        
    case {'DelayRemaining' 'PercentDelayRemaining' 'TimeToQuit' 'TimeToLeaveFood'}
        nLaps = length(LapStart);
        LapTime = vt.x;
        LapTime.D(:) = nan;
        timeData = LapTime.range;
        for iL = 1:nLaps 
            if isnan(LapStart(iL)) % Delay fields will have nan for skip laps
                continue
            end
             validIdx = timeData >= LapStart(iL) & timeData <= LapEnd(iL);
            tempTime = timeData(validIdx);
            currTime = LapData(iL) - (tempTime - LapStart(iL));
            if strcmp(Variable,'PercentDelayRemaining')
                currTime = currTime/LapData(iL)*100;
            end
            LapTime.D(validIdx) = currTime;
        end        
        outputTSD = LapTime;
        switch Variable
            case 'PercentDelayRemaining'
                binEdges = [0 100];
            case 'DelayRemaining'
                binEdges = [0 30];
            case 'TimeToQuit'
                binEdges = [0 25];
            case 'TimeToLeaveFood'
                binEdges = [0 34];
        end
        
           
        
        % These are made into TSDs from existing TSDs
    case 'Random'
        randMax = 1000;
        RandInput = randi(randMax,length(vt.x.range),1);
        outputTSD = tsd(vt.x.range,RandInput);
        binEdges = [0 randMax];
    case 'RandomSmooth'
        randMax = 1000;
        smoothSigma = .1; %Sec
        smoothWindow = 1;
        RandInput = randi(randMax,length(vt.x.range),1);
        temp = tsd(vt.x.range,RandInput);
        outputTSD = smooth(temp, smoothSigma, smoothWindow);
        binEdges = [min(outputTSD.D) max(outputTSD.D)];
        
    case 'SessionTime'
        outputTSD = tsd(vt.x.range,vt.x.range-vt.x.starttime);
        binEdges = [0 3600];
        
         
        
        %%% These are made into TSDs from raw info %%%
    case 'RandomLap' % Random value for each lap but constant for the given lap
        RandInput = LapData;
        outputTSD= BuildEventTSD(LapStart,RandInput,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 max(RandInput(:))];
        
    case 'IdPhi'
        IdPhi = log(LapData); % Pass in raw IdPhi
        outputTSD= BuildEventTSD(LapStart,IdPhi,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 9];
    case 'zIdPhi'
        IdPhi = LapData; % Pass in raw IdPhi
        zIdPhi = (IdPhi - nanmean(IdPhi(:)))./nanstd(IdPhi(:));
        outputTSD= BuildEventTSD(LapStart,zIdPhi,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [-2 10];
    case 'AvgDPhi'
        AvgDPhi = LapData;
        outputTSD = BuildEventTSD(LapStart,AvgDPhi,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 8];
        
    case 'Delay'
        Delay = LapData;
        outputTSD= BuildEventTSD(LapStart,Delay,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [.5 30.5];
    case 'Value'
        Value = LapData;% Pass in fully computed value
        outputTSD= BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [-25 25];  
    case 'OfferDifficulty' % 
        Difficulty = -abs(LapData);% Pass in fully computed value
        outputTSD= BuildEventTSD(LapStart,Difficulty,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [-25 0];
    case 'OfferEase'
        Ease = abs(LapData);% Pass in fully computed value
        outputTSD= BuildEventTSD(LapStart,Ease,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 25];
        
    case {'TotalEarns' 'EarnedRewards'}
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 160];
    case {'PercentTotalEarns'}
        Value = LapData/max(LapData)*100;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));        
        binEdges = [0 100];
    case 'LapsFromLastEarn'
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [.5 10.5];
        
    case 'ReactionTime'
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 7.5];
    case {'ConsumptionTime' 'LingerTime'}
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [4 34];
    case 'PrctLingerTime'
         Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [3 30];        
    case 'QuitTime'
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [0 25];
        
        %%% Categorical Data %%%
    case {'SiteNumber' 'NextSite' 'Rank' 'SiteRank' 'NextRank' 'PreviousSite' 'PreviousRank'}
        Value = LapData; % Location for site, Rank for rank
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [.5 4.5];
        nBins = 4;
        
    case {'Skip' 'Accept' 'Earn' 'Quit' 'NoEarnLap' 'EarnLap' 'GoodDeal' 'BadDeal' 'PreviousAccept' 'NextAccept' 'EconomicChoice'}
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [-.5 1.5];
        nBins = 2;
        
    case 'VTE'
        IdPhi = LapData;
        [th] = HalfGaussVteThresh(IdPhi);
        VTELaps = IdPhi > th;
        outputTSD= BuildEventTSD(LapStart,VTELaps,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [-.5 1.5];
        nBins = 2;
    case 'VTE_BySite' % Threshold calculated for each site independently
        IdPhi = LapData;
        VTELaps = nan(size(IdPhi));
        
        zoneIDs = unique(LapSite);
        for iZ = 1:length(zoneIDs)
            LapsToUse = LapSite == zoneIDs(iZ);
            [th] = HalfGaussVteThresh(IdPhi(LapsToUse));
            VTELaps(LapsToUse) = IdPhi(LapsToUse) > th;
        end
        outputTSD= BuildEventTSD(LapStart,VTELaps,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [-.5 1.5];
        nBins = 2;
        
    case 'OZRunSpeed'
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [5 50];
        
    case {'Choice' 'Decision'}
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [.5 3.5];
        nBins = 3;
        
    case 'MazeChunk'
        Value = LapData;
        outputTSD = BuildEventTSD(LapStart,Value,'fullTS',fullTS,'endTime',max(LapEnd(:)));
        binEdges = [.5 4.5];
        nBins = 4;
        
    otherwise
        error(sprintf('Case (%s) not implimented yet',Variable))
end
