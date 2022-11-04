function [LapStart, LapEnd, LapSite, LapAdvance, LapPreStart,...
    OZEnter,OZExit, WZEnter,WZExit, LZEnter,LZExit, AZEnter,AZExit, FullEventTimes, WZExit_Adj, LZExit_Adj, AZEnter_Adj] ...
    = identifyRRowLapTimes(LapData,OfferZones,WaitZones,varargin)

% Takes a combined behavioral output stucture generated from
% 'collectRRowDataToAnalyze' and extracts out a series of relevant task
% time points (StartTimes, EndTimes of various task points).
%
%
% INPUT
% LapData - An aggregate structure generated from
%       'collectRRowDataToAnalyze' that lists most/all RRow behavioral/task
%       results. Each entry is in a Sess x Lap x Site format. In the event
%       that this is only a single session it can be passed in as Lap x
%       Site and the code manages it approprietly.
% OfferZones - What site number are the offer zones in the FullLapData.
%       Standard conditions (4 site RRow) these would be 5:8
% WaitZones - What site number are the offer zones. Standard would be 1:4.
%
% epsi - Small time step to pullback exit times from entry times such that
%       they do not overlap. This should be >> matlab resolution and < any
%       cammera/phys data sampling
% sessNum - Can designate only a subset of sessions (most usefully 1) to
%       get output for.
% squeezeD1 - In the event that only a sinlge session is extracted you can
%       squeeze out the 1st dimension (session). Default true.
%
% OUTPUTS
% Series of Sess x Lap x Site matrices of consolidated lap times.

% GWD 2019


epsi = 1e-10;
tsdStep = 1e-3;
sessNum = [];
squeezeD1 = false;
process_varargin(varargin);


nDim =  ndims(LapData.Decision);
if nDim == 2
    fieldEntries = fieldnames(LapData);
    squeezeD1 = true;
    for iF = 1:length(fieldEntries)
        LapData.(fieldEntries{iF}) = permute(LapData.(fieldEntries{iF}),[3 1 2]);
    end    
end
[nSess, nLaps, nZones] = size(LapData.Decision);
if isempty(sessNum)
    sessNum = 1:nSess;
end


% Entering the offer zone
LapStart = LapData.EnteringZoneTime(sessNum,:,OfferZones);

% Instant before the next offer zone entry (dont want next time to be kept
% in the restriction). Also step back 1.5 tsd steps so still pulls from the
% same entry in lapwise based TSDs
LapEnd = advanceRRowData(LapStart,-1) - epsi - 1.5*tsdStep;
for iS = 1:length(sessNum)
    LapEnd(sessNum(iS),end) = nan; % Kill the start of the first lap that gets looped to the back
end

% Site/Restaurant number
LapSite = LapData.ZoneNum(sessNum,:,OfferZones) - nZones;

% Point at which everything is done and advancing to the next RR
AdvanceAccept = LapData.ExitZoneTime(sessNum,:,WaitZones); % Exits from Accepted Waits (Earn/Quit)
AdvanceSkip = LapData.ExitZoneTime(sessNum,:,OfferZones); % Exits from OfferZones
AdvanceSkip(LapData.AcceptOffer(sessNum,:,OfferZones)==1) = nan; % Get rid of OfferZone exits when accepting
LapAdvance = nanmean(cat(4,AdvanceAccept,AdvanceSkip),4); %Nanmean to take only the real value

% Point in the preceeding lap at which everything is done (in the advance
% zone)
LapPreStart = advanceRRowData(LapAdvance,1) + epsi;


% Simple set of enter/exit times for each of the 4 maze chunks. In several
% cases these are redundent with some of the upper values but will make
% life easier in some cases to have it clean/clear.
OZEnter = LapStart;
OZExit = LapData.ExitZoneTime(sessNum,:,OfferZones) - epsi;

WZEnter = LapData.EnteringZoneTime(sessNum,:,WaitZones);
earnTimes = LapData.FeederTimes(sessNum,:,WaitZones);
quitTimes = LapData.ExitZoneTime(sessNum,:,WaitZones);
quitTimes(LapData.EarnOffer(sessNum,:,WaitZones)==1) = nan;
% Nanmean to take only the valid entries while maintaining shape
WZExit = nanmean(cat(4,earnTimes,quitTimes),4) - epsi;


LZEnter = LapData.FeederTimes(sessNum,:,WaitZones);
lingerLeaves = LapData.ExitZoneTime(sessNum,:,WaitZones);
lingerLeaves(LapData.QuitOffer(sessNum,:,WaitZones)==1) = nan;
LZExit = lingerLeaves - epsi;


AZEnter = LapAdvance;
AZExit = LapEnd;

if isfield(LapData,'AdjWZExitTime')
    % Adjust the respective Exit (WZ/LZ) and Entry (AZ) for crossing the OZ
    % threshold location. Note that WZExit is only adapted for Quit decisions
    % (Earns will have WZExit flowing into LZEntry). LZExit is adapted for
    % quits.
    AdjWZExit = LapData.AdjWZExitTime(sessNum,:,WaitZones);
    
    % Take the latest AZEnter times. For skips these will be concurrent,
    % for Accepts (earns or quits) the Adj will be after the standard (WZ
    % zone exit vs OZ zone exit)
    AZEnter_Adj = max(cat(4,AZEnter,AdjWZExit),[],4);
    
    % For Earned offers take the adjusted WZExit zone time and remove any
    % events without LZExit times (Quits/Skips)
    LZExit_Adj = max(cat(4,LZExit,AdjWZExit - epsi),[],4);
    LZExit_Adj(isnan(LZExit)) = nan;
    
    % WZExit for Earns is the same (those adjust the LZExit), but quits we
    % are going to take the Adj zone exit time
    WZExit_Adj = WZExit;
    WZExit_Adj(LapData.QuitOffer(sessNum,:,WaitZones)==1) = AdjWZExit(LapData.QuitOffer(sessNum,:,WaitZones)==1) - epsi;
else
    AZEnter_Adj = nan(size(WZExit));
    LZExit_Adj = nan(size(WZExit));
    WZExit_Adj = nan(size(WZExit));
end

% Boil things down to a single/subset of sessions if desired
if length(sessNum) == 1 && squeezeD1
    LapStart = reshape(LapStart,nLaps,nZones);
    LapEnd = reshape(LapEnd,nLaps,nZones);
    LapSite = reshape(LapSite,nLaps,nZones);
    LapAdvance = reshape(LapAdvance,nLaps,nZones);
    LapPreStart = reshape(LapPreStart,nLaps,nZones);
    
    OZEnter = reshape(OZEnter,nLaps,nZones);
    OZExit = reshape(OZExit,nLaps,nZones);
    WZEnter = reshape(WZEnter,nLaps,nZones);
    WZExit = reshape(WZExit,nLaps,nZones);
    
    LZEnter = reshape(LZEnter,nLaps,nZones);
    LZExit = reshape(LZExit,nLaps,nZones);
    AZEnter = reshape(AZEnter,nLaps,nZones);
    AZExit = reshape(AZExit,nLaps,nZones);
    
    AZEnter_Adj = reshape(AZEnter_Adj,nLaps,nZones);
    LZExit_Adj = reshape(LZExit_Adj,nLaps,nZones);
    WZExit_Adj = reshape(WZExit_Adj,nLaps,nZones);
end

% Single structure with all fields so that it can be dynamically referenced
FullEventTimes = struct;
FullEventTimes.LapStart = LapStart;
FullEventTimes.LapEnd = LapEnd;
FullEventTimes.LapSite = LapSite;
FullEventTimes.LapAdvance = LapAdvance;
FullEventTimes.LapPreStart = LapPreStart;

FullEventTimes.OZEnter = OZEnter;
FullEventTimes.OZExit = OZExit;
FullEventTimes.WZEnter = WZEnter;
FullEventTimes.WZExit = WZExit;

FullEventTimes.LZEnter = LZEnter;
FullEventTimes.LZExit = LZExit;
FullEventTimes.AZEnter = AZEnter;
FullEventTimes.AZExit = AZExit;

FullEventTimes.AZEnter_Adj = AZEnter_Adj;
FullEventTimes.LZExit_Adj = LZExit_Adj;
FullEventTimes.WZExit_Adj = WZExit_Adj;