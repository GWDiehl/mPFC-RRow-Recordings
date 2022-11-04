
function [LapData, SessionData, IdPhi] = collectRRowDataToAnalyze(varargin)

% Collects up RRow behaivoral data and stores it in a single strucutre with
% each field having Sess x Laps x Sites matrix of entries.
% Loads behavioiral data for each lap (LapData), relevent data for each
% session (SessionData), and IdPhi data for each lap (IdPhi).
%
% In general, the FULL slate of data has been processed and saved for
% future recall as it can take a while to run through (~1hr). Subsequently
% this full data can be loaded and segmented down to just what is wanted.
% Segmentation is done based on matching fields to dataDef file.
%
% If saved data does not exist it will be processed de-novo from the full
% dataDef info.
%
%
% Optional Inputs
% OverwriteData - Logical that will recalculate and overwrite all data as
%       opposed to loading anything that is saved.
% Task - Str of which task to process. (Currently only 'RRow' supported).
%       Different tasks should be procssed and saved into different files
%       as they presumibly have different conditions/needs/considerations
%
% Restrictions - Matching cell array of field and values used to do a
%   logical AND restriction of full data set.
% restrictionField - Cell array of what field in dataDef you would like to
%       use to restrict the full data structure.
% restrictionValue - Cell array of what value/entry/str/etc. you would like
%       to fall within/match for restriction. Can be a single str, cell
%       array of strs, single value to match, or pair of values to fall
%       within.

% GWD 2019

[RootDir, DataDefName, BehavDir, TaskSuffix, PathSuffix,...
    ~, ~, ~, ~, ~, PromotionDir] = ID_GWD_Dirs;


SessDataSuffix = '';
LapDataSuffix = '';
IdPhiDataSuffix = '';

OverwriteData = 0;
Task = 'RRow'; % or 'RRow4x20' - Not yet implemented

restrictionField = cell(0);
restrictionValue = cell(0);

% Threshold for identifying "running" to be used in calculaing avg run speed to stationary time.
runSpeedThreshold = 5; 

process_varargin(varargin);



%% Load datadef and file names

dataDef = load([RootDir,DataDefName]);


processedSessions = cellfun(@(x,y) [RootDir,BehavDir,x,'\',y,'\'],dataDef.RatID,dataDef.Directory,'uniformoutput',0);
filePrefix = dataDef.SSN;
fullTaskFN = cellfun(@(x) [x,TaskSuffix],filePrefix,'UniformOutput',0);
fullPathFN = cellfun(@(x) [x,PathSuffix],filePrefix,'UniformOutput',0);


%% Collect up the behavioral data or load from save

correctTask = ismember(dataDef.Task,Task);
subDataDef = subsetFromDataDef(dataDef,correctTask);

% If the user did not input any prefered suffix default based on Task
if isempty([SessDataSuffix LapDataSuffix IdPhiDataSuffix])
    switch Task
        case 'RRow'
            % All good, no suffix
        case 'RRow4x20'
            SessDataSuffix = '_4x20';
            LapDataSuffix = '_4x20';
            IdPhiDataSuffix = '_4x20';
    end
end

% Task Data
if exist([RootDir, BehavDir 'SessionData_RRow',SessDataSuffix,'.mat'],'file') && ~OverwriteData
    SessionData = load([RootDir, BehavDir 'SessionData_RRow',SessDataSuffix]);
    LapData = load([RootDir, BehavDir, 'LapData_Behav_RRow',LapDataSuffix]);    
else
    switch Task
        case 'RRow'
            [SessionData, LapData] = collectRRowBehaviorData(subDataDef,...
                processedSessions(correctTask),fullPathFN(correctTask),fullTaskFN(correctTask),'runSpeedThreshold',runSpeedThreshold,'trackingDataSwap',{BehavDir, PromotionDir});    
        case 'RRow4x20'
            error('Write in this case')
    end    
    save([RootDir, BehavDir 'SessionData_RRow',SessDataSuffix],'-struct','SessionData')
    save([RootDir, BehavDir, 'LapData_Behav_RRow',LapDataSuffix],'-struct','LapData')
end

%IdPhi Data
if exist([RootDir,BehavDir,'IdPhi_RRow',IdPhiDataSuffix,'.mat'],'file') && ~OverwriteData
    IdPhi = load([RootDir,BehavDir,'IdPhi_RRow',IdPhiDataSuffix]);
else
    switch Task
        case 'RRow'
            [~, IdPhi] = collectRRowIdPhi(subDataDef,processedSessions(correctTask),fullPathFN(correctTask),fullTaskFN(correctTask));
        case 'RRow4x20'
            error('Write in this case')
    end
    save([RootDir,BehavDir,'IdPhi_RRow',IdPhiDataSuffix],'-struct','IdPhi')
end


% If we loaded the data, make sure that we segment the dataDef so that it
% matches

assert(isequal(size(LapData.EnteringZoneTime,1),size(SessionData.SessDate,1),size(IdPhi.IdPhi,1)),'Your behaviors do not all match, you have a problem')
validDefSess = cellfun(@(x,y) ismember(x,SessionData.SessDate) & ismember(y,SessionData.RatID),subDataDef.Date,subDataDef.RatID);
subDataDef = subsetFromDataDef(subDataDef,validDefSess);

assert(isequal(SessionData.SessDate(:),subDataDef.Date(:)) && isequal(SessionData.RatID(:),subDataDef.RatID(:)),'Your Session listings do not match, you have a problem.')


%% Restrict Sessions to pass out based on dataDef restrictions

sessionsToUse = true(1,length(subDataDef.RatID));
for R = 1:length(restrictionField)
    % Restricting based on a str or set of strs
    if ischar(restrictionValue{R}) || iscell(restrictionValue{R})
        sessionsToUse = sessionsToUse & ismember(subDataDef.(restrictionField{R}),restrictionValue{R});
        
        % Restricting based on matching a single value
    elseif isnumeric(restrictionValue{R}) && legnth(restrictionValue{R}) == 1
        sessionsToUse = sessionsToUse & subDataDef.(restrictionField{R}) == restrictionValue{R};
        
        % Restricting based on being within a pair of values
    elseif isnumeric(restrictionValue{R}) && legnth(restrictionValue{R}) == 2
        sessionsToUse = sessionsToUse & iswithin(subDataDef.(restrictionField{R}),restrictionValue{R}(1),restrictionValue{R}(2));        
    end
end

if any(~sessionsToUse)
    LapData = subsetStructByIndx(LapData,sessionsToUse);
    SessionData = subsetStructByIndx(SessionData,sessionsToUse);
    IdPhi = subsetStructByIndx(IdPhi,sessionsToUse);
end


%% Add in anything else you may want but does not come directly from the origianl data (thus may be modified/tweaked or such)

[nSess, nLaps, nZones] = size(LapData.ZoneDelay);

LapData.Value_H = repmat(reshape(SessionData.Threshold.Global.Heaviside,nSess,1,nZones/2),1,nLaps,2) - LapData.ZoneDelay;
LapData.Value_S = repmat(reshape(SessionData.Threshold.Global.Sigmoid,nSess,1,nZones/2),1,nLaps,2) - LapData.ZoneDelay;
LapData.TotalSiteTime = nansum(cat(4,LapData.ZoneTime(:,:,1:(nZones/2)),LapData.ZoneTime(:,:,(nZones/2)+1:end)),4);
LapData.TotalSiteTime(isnan(LapData.ZoneDelay(:,:,(nZones/2)+1:end))) = nan;

thresholdVals = SessionData.Threshold.Global.Heaviside;
nRewardsEarned = reshape(nansum(LapData.EarnOffer(:,:,1:nZones/2),2),nSess,nZones/2);
[FeederRank,RankOrder] = ExtractRRowRankOrdering(cat(3,thresholdVals,nRewardsEarned));
LapData.SiteRank = permute(repmat(FeederRank,1,2,nLaps),[1 3 2]);


IdPhiFields = fieldnames(IdPhi);
for iF = 1:length(IdPhiFields)
    LapData.(IdPhiFields{iF}) = IdPhi.(IdPhiFields{iF});
end

sessStartTime = arrayfun(@(x) min(LapData.EnteringZoneTime(x,:)),1:nSess)';
sessStartTime = repmat(sessStartTime,1,nLaps,nZones);
LapData.ElapsedSessTime = LapData.EnteringZoneTime - sessStartTime;



