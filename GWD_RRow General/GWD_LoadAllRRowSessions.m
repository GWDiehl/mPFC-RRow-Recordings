function [SSNs,vt,R,S,S_FNs] = GWD_LoadAllRRowSessions(varargin)


% Loads ALL RRow data (SSN, tracking, Behavior R, Spikes, SpikeFN) for GWD
% mPFC project. Data can be restricted based on subselections of dataDef
% file using TaskOI and ExpOI. Default takes all sessions with Phys data
% ('Recording-Stable') in the standard RRow task (RRow).
%
% Select data sets can be gotten by restricting with TaskOI, ExpOI, or
% SSNsToUse
%
% Data output as a cell array with one entry for each session
% (corresponding across data).

% GWD 2019

TaskOI = 'RRow';
ExpOI = 'Recording-Stable';
dataToLoad = {'Task' 'Path' 'Spikes'};
SSNsToUse = [];
restrictToTaskTime = 1;
process_varargin(varargin);

%%

% Directories
[rootDir, DataDefName] = ID_GWD_Dirs;

% Data Def loading

dataDef = load([rootDir,DataDefName]);

if isempty(SSNsToUse)
    DataDefOI = ismember(dataDef.Task,TaskOI) & ismember(dataDef.ExpType,ExpOI);
else
    DataDefOI = ismember(dataDef.SSN,SSNsToUse);
end

% Grab the sub directories
subDataDef = subsetFromDataDef(dataDef,DataDefOI);

SSNs = subDataDef.SSN';
nSess = length(SSNs);

vt = cell(nSess,1);
R = cell(nSess,1);
S = cell(nSess,1);
S_FNs = cell(nSess,1);

for iS = 1:nSess
    
    [vt{iS}, R{iS}, S{iS}, S_FNs{iS}] = GWD_LoadRRowSession(SSNs{iS},dataToLoad,'restrictToTaskTime',restrictToTaskTime);
    
end

