function [vt, R, S, S_FNs,S_FNs_Ext,LFP,LFP_FNs,LFP_FNs_Ext] = GWD_LoadRRowSession(SSN,dataToLoad,varargin)

% General loading engine for GWD RRow data set. Will take a single session
% input SSN as well as a cell array of str indicting which data you would
% like loaded {'Path' 'Task' 'Spikes' 'SpikeNames'}
%
% Tacking data (vt) and spikes data (S) are returned as standard ts/tsd
% objected. Task data (R) is a cleaned up version of the original output
% from Matlab
%
% Optional inputput 'restrictToTaskTime' will restrict tracking/spiking
% data to time on track as defined by keys file.

% GWD 2019

% NOTE: LFP Loading not currently supported

%%

if iscell(SSN)
    nSess = length(SSN);
    if nSess == 1
        SSN = SSN{:};
    else
        error('Loop through sessions on the outside, just pass in a single SSN');
    end
end

[RootDir, ~, BehavDir, TaskFileRoot, PathFileRoot, KeysFileRoot, SpikesDir, LFPDir, ChannelMapDir] = ID_GWD_Dirs;

TFileRoot = '.t';
LFPDataRoot = '.lfp';
LFPTimeRoot = '-LfpTime.dat';

restrictToTaskTime = 1;

LFPsToLoad = {};


process_varargin(varargin);

if ~iscell(LFPsToLoad)
    LFPsToLoad = {LFPsToLoad};
end

% Ensure that it is a single vector
LFPsToLoad = LFPsToLoad(:);


%% Do stuff

RatID = SSN(1:4);

TaskFN = fullfile(RootDir,BehavDir,RatID,SSN,[SSN TaskFileRoot]);
PathFN = fullfile(RootDir,BehavDir,RatID,SSN,[SSN PathFileRoot]);
keysFN = [strrep(SSN,'-','_') KeysFileRoot];

vt = [];
R = [];
S = cell(0);
S_FNs = cell(0);
S_FNs_Ext = cell(0);
LFP = cell(0);
LFP_FNs = cell(0);
LFP_FNs_Ext = cell(0);

if restrictToTaskTime
    origDir = pwd;
    cd(fileparts(PathFN))
    [~,keysFunc] = fileparts(keysFN);
    keys = evalKeysFile(keysFunc);
    startTask = keys.TimeOnTrack;
    endTask = keys.TimeOffTrack;
    cd(origDir);
end

if ismember('Task',dataToLoad)
    R = load(TaskFN);
    if isfield(R,'R')
        R = R.R;
    end
    
end
if ismember('Path',dataToLoad)
    vt = load(PathFN);
    if isfield(vt,'vt')
        vt = vt.vt;
    end
    if restrictToTaskTime
        pathFeatures = fields(vt);
        for iF = 1:length(pathFeatures)
            vt.(pathFeatures{iF}) = vt.(pathFeatures{iF}).restrict(startTask,endTask);
        end
    end
    
end

if ismember('Spikes',dataToLoad) || ismember('SpikeNames',dataToLoad)
    
    SpikesDir = fullfile(RootDir,SpikesDir,RatID,SSN); 
    fns = FindFiles(['*',TFileRoot],'StartingDirectory',SpikesDir,'CheckSubdirs',0);
    [~, S_FNs,S_FNs_Ext] = cellfun(@(x) fileparts(x),fns,'UniformOutput',0);
    
    if ismember('Spikes',dataToLoad)
        S = LoadSpikes(fns);
        if restrictToTaskTime
            S = cellfun(@(x) x.restrict(startTask,endTask),S,'uniformoutput',0);
        end
    end
end

if ismember('LFP',dataToLoad) || ismember('LFPNames',dataToLoad)
    LFPDir = fullfile(RootDir,LFPDir,RatID,SSN);
    if isempty(LFPsToLoad)
        LFPsToLoad = FindFiles(['*',LFPDataRoot],'StartingDirectory',LFPDir,'CheckSubdirs',0);
    else
        if ~contains(LFPsToLoad{1},LFPDir)
            if contains(LFPsToLoad{1},SSN)
                LFPsToLoad = cellfun(@(x) [LFPDir,'\',x,LFPDataRoot],LFPsToLoad,'UniformOutput',0);
            else
                LFPsToLoad = cellfun(@(x) [LFPDir,'\',SSN,'-',x,LFPDataRoot],LFPsToLoad,'UniformOutput',0);
            end
        end
    end
    
    [~, LFP_FNs,LFP_FNs_Ext] = cellfun(@(x) fileparts(x),LFPsToLoad,'UniformOutput',0);    
    
    if ismember('LFP',dataToLoad)
        time_FN = FindFiles(['*',LFPTimeRoot],'StartingDirectory',LFPDir,'CheckSubdirs',0);
        assert(length(time_FN)==1,'You have more than one LFP time file in the directory. That is a problem')
        
        timeData = LoadIntanLfpTime(time_FN{1});
        
        LFP = cell(length(LFPsToLoad),1);
        for iN = 1:length(LFPsToLoad)
            lfpData = LoadIntanLfp(LFPsToLoad{iN});
            
            LFP{iN} = tsd(timeData,lfpData);
            
            if restrictToTaskTime
                LFP{iN} = restrict(LFP{iN},startTask,endTask);
            end
        end
    end
    
end

