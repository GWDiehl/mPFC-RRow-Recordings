
function [SessionData, LapData] = collectRRowIdPhi(DataDef,SessionDirs,PathFileName,BehaviorFileName,varargin)


maxLaps = 125;

numSites = 4;
WaitZones = 1:numSites;
OfferZones = numSites+1:2*numSites;
sitesOfInterest = [WaitZones OfferZones];

% Amount of time (sec) to use from OZ entry or WZ exit in IdPhi calc
% Alternativly use inf to just work based on opposing zone entry/exit time.
OZIdPhiTime = inf;
WZIdPhiTime = 3;

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

for s = 1:length(SessionDirs)
    
    vt = load([SessionDirs{s},'\',PathFileName{s}]);
    R = load([SessionDirs{s},'\',BehaviorFileName{s}]);
    if isfield(R,'R')
        R = R.R;
    end
    if isfield(vt,'vt')
        vt = vt.vt;
    end
    
    groupedData = groupRRowBySite(R, 'sitesOfInterest', sitesOfInterest, 'maxLaps', maxLaps, 'fieldsOfInterest', {'EnteringZoneTime' 'ExitZoneTime'});
    
    EnterTimes = groupedData.EnteringZoneTime;
    ExitTimes = groupedData.ExitZoneTime;
    
    ExitTimes(:,OfferZones) = min(ExitTimes(:,OfferZones),EnterTimes(:,OfferZones)+OZIdPhiTime);
    ExitTimes(:,WaitZones) = min(ExitTimes(:,WaitZones),EnterTimes(:,WaitZones)+WZIdPhiTime);
    
    [IdPhi, AvgDPhi] = compute_IdPhi_AvgdPhi(vt.x,vt.y,EnterTimes(:),ExitTimes(:));
    
    IdPhi = reshape(IdPhi,size(EnterTimes));
    AvgDPhi = reshape(AvgDPhi,size(EnterTimes));
    
    LapData.IdPhi(s,:,:) = IdPhi;
    LapData.AvgDPhi(s,:,:) = AvgDPhi;
    
    fprintf('Processed IdPhi for %s-%s \n',SessionData.RatID{s},SessionData.SessDate{s})
end


