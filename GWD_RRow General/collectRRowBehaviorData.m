
function [SessionData, LapData] = collectRRowBehaviorData(DataDef,SessionDirs,PathFileName,BehaviorFileName,varargin)


maxLaps = 125;
numDays = length(SessionDirs);



numSites = 4;
WaitZones = 1:numSites;
OfferZones = numSites+1:2*numSites;
sitesOfInterest = [WaitZones OfferZones];

% This is by far the longest calculation so if not needed just don't do it
calcRunSpeed = 1;
runSpeedThreshold = 5;

% How much to we want to goose the WZ exit time forward when finding the
% time that the rat pass out of the OZ (done to make sure that the start
% point of the path is really IN the OZ).
adjWZExitTime = .15;

trackingDataSwap = {'',''};
process_varargin(varargin);


%% Initialize Groupings


SessionData = struct;
LapData = struct;

SessionData.SessDate = DataDef.Date';
SessionData.RatID = DataDef.RatID';
SessionData.ExpType = DataDef.ExpType';

SessionData.TaskPhase = DataDef.TaskPhase';
SessionData.PhaseDay = DataDef.PhaseDay';
SessionData.TrainingDay = DataDef.TrainingDay';


SessionData.Threshold.Global.Heaviside = nan(numDays,numSites);
SessionData.Threshold.OfferZone.Heaviside = nan(numDays,numSites);
SessionData.Threshold.WaitZone.Heaviside = nan(numDays,numSites);

SessionData.Threshold.Global.Sigmoid = nan(numDays,numSites);
SessionData.Threshold.OfferZone.Sigmoid = nan(numDays,numSites);
SessionData.Threshold.WaitZone.Sigmoid = nan(numDays,numSites);

SessionData.Threshold.Global.Sigmoid_Slope = nan(numDays,numSites);
SessionData.Threshold.OfferZone.Sigmoid_Slope = nan(numDays,numSites);
SessionData.Threshold.WaitZone.Sigmoid_Slope = nan(numDays,numSites);


for s = 1:length(SessionDirs)
    
   
    R = load([SessionDirs{s},'\',BehaviorFileName{s}]);
    if isfield(R,'R')
        R = R.R;
    end
     
    useOfferZone = SessionData.TaskPhase(s) > 2; % Task phase of 3+ has offer zone
    
    groupedData = groupRRowBySite(R, 'sitesOfInterest', sitesOfInterest, 'maxLaps', maxLaps);
    
    %% Add in extra things
    
    % Add Accept/Earn Fields
    groupedData.AcceptOffer = nan(maxLaps,length(sitesOfInterest));
    groupedData.AcceptOffer(groupedData.SkipOffer == 0) = 1;
    groupedData.AcceptOffer(groupedData.SkipOffer == 1) = 0;
    groupedData.EarnOffer = nan(maxLaps,length(sitesOfInterest));
    groupedData.EarnOffer(groupedData.QuitOffer == 0) = 1;
    groupedData.EarnOffer(groupedData.QuitOffer == 1) = 0;
    
    % Add Zone info
    groupedData.ZoneNum = repmat(sitesOfInterest,maxLaps,1);
    groupedData.ZoneID = cell(maxLaps,length(sitesOfInterest));
    groupedData.ZoneID(1:maxLaps,WaitZones) = {'WaitZone'};
    groupedData.ZoneID(1:maxLaps,OfferZones) = {'OfferZone'};
    
    % Add decision info
    groupedData.Decision = cell(maxLaps,numSites);
    groupedData.Decision(:) = {'N/A'};
    groupedData.Decision(groupedData.SkipOffer(:,OfferZones) == 1) = {'Skip'};
    groupedData.Decision(groupedData.QuitOffer(:,WaitZones) == 1) = {'Quit'};
    groupedData.Decision(groupedData.EarnOffer(:,WaitZones) == 1) = {'Earn'};
    
    
    groupedData.ZoneTime = groupedData.ExitZoneTime - groupedData.EnteringZoneTime;
    
    
    %% Thresholds for each reward site (done independently)
    
    for z = 1:numSites
        % Global Threshold: Was there a reward given for each offer?
        
        % Take the delay from Wait or Offer Zone. Always identical but
        % could be nan for skips/quit/etc.
        Delay = groupedData.ZoneDelay(:,OfferZones(z));
        % Was the offer earned, regardless of if it was a skip (nan in wait zone)
        Rewarded = groupedData.EarnOffer(:,WaitZones(z)) == 1;
        
        SessionData.Threshold.Global.Heaviside(s,z) = RRheaviside(Delay,Rewarded);
        [Threshold,Slope] = RRowSigmoidFit(Delay,Rewarded);
        SessionData.Threshold.Global.Sigmoid(s,z) = Threshold;
        SessionData.Threshold.Global.Sigmoid_Slope(s,z) = Slope;
        
        
        % WaitZone Threshold: Was reward earned given the offers the were
        % intially accepted (Ignore Skips)?
        Delay = groupedData.ZoneDelay(:,WaitZones(z)); % Info from Wait Zones
        Earn = groupedData.EarnOffer(:,WaitZones(z));
        
        SessionData.Threshold.WaitZone.Heaviside(s,z) = RRheaviside(Delay,Earn);
        [Threshold,Slope] = RRowSigmoidFit(Delay,Earn);
        SessionData.Threshold.WaitZone.Sigmoid(s,z) = Threshold;
        SessionData.Threshold.WaitZone.Sigmoid_Slope(s,z) = Slope;
        
        
        if useOfferZone
            % OfferZone Threshold: Was the offer accepted for each offer?
            Delay = groupedData.ZoneDelay(:,OfferZones(z));
            Accept = groupedData.AcceptOffer(:,OfferZones(z)); % Info from Offer Zones
            
            SessionData.Threshold.OfferZone.Heaviside(s,z) = RRheaviside(Delay,Accept);
            [Threshold,Slope] = RRowSigmoidFit(Delay,Accept);
            SessionData.Threshold.OfferZone.Sigmoid(s,z) = Threshold;
            SessionData.Threshold.OfferZone.Sigmoid_Slope(s,z) = Slope;
        end
        
    end
    
    %% Extract out the Adjusted WZTimes
    
    vt = load([strrep(SessionDirs{s},trackingDataSwap{1},trackingDataSwap{2}),'\',PathFileName{s}]);
    if isfield(vt,'vt')
        vt = vt.vt;
    end
    
    groupedData.AdjWZExitTime = nan(maxLaps,length(sitesOfInterest));
    
    taskZones = R.params.world.Zones;
    % The original run code had the OZ as 1:4 and WZ as 5:8; GWD Revamp
    % made it as OZ 5:8 and WZ as 1:4 to prioratize being in the WZ
    if strcmp(R.params.RunCode,'RRow_2zone')
        assert(length(taskZones) == 8, 'Why do you not have 8 zones in this run config??')
        tempZones = taskZones;
        taskZones(OfferZones) = tempZones(WaitZones);
        taskZones(WaitZones) = tempZones(OfferZones);
    end
    
    [~,~,~,~,~,~,~,~,WZExit,~,LZExit,~,~,~]...
        = identifyRRowLapTimes(groupedData,OfferZones,WaitZones);
    
    selectedTime = max(cat(3,WZExit,LZExit),[],3)' + adjWZExitTime;
    OZExitTimes = nan(size(selectedTime));
    validTimes = ~isnan(selectedTime);
    currZone = repmat([1:numSites],maxLaps,1)';
    
    % Make sure that the tracking is all positive values
    tempVT = vt;
    tempVT.y.D = abs(tempVT.y.D);
    tempVT.x.D = abs(tempVT.x.D);
    
    temp = extractEarnExits(selectedTime(validTimes),tempVT.x,tempVT.y,taskZones,currZone(validTimes),'ZoneSpread',8,'convertNlx',0);
    OZExitTimes(validTimes) = temp;
    
    groupedData.AdjWZExitTime(:,WaitZones) = OZExitTimes';
    
    %% Run Speed
    
    if calcRunSpeed
        vt = load([SessionDirs{s},'\',PathFileName{s}]);
        if isfield(vt,'vt')
            vt = vt.vt;
        end
        sampleRate = 1/median(diff(vt.x.range));
        
        RunSpeed = arrayfun(@(x,y) vt.v.restrict(x,y).D,groupedData.EnteringZoneTime,groupedData.ExitZoneTime,'uniformoutput',0);
        groupedData.RunSpeed = cellfun(@(x) nanmedian(x(x>runSpeedThreshold)),RunSpeed);
        groupedData.PauseTime = cellfun(@(x) sum(x<runSpeedThreshold)/sampleRate,RunSpeed);
        groupedData.PauseTime(isnan(groupedData.ZoneTime)) = nan;
    end
    
    %% Colelct up everything into a big structure across sessions
    entries = fields(groupedData);
    for f = 1:length(entries)
        LapData.(entries{f})(s,:,:) = groupedData.(entries{f});
    end
    
    fprintf('Processed Behavior for %s-%s \n',SessionData.RatID{s},SessionData.SessDate{s})
end



