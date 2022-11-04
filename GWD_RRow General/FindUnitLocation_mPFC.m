function [stdLocation, netDepth, peakChan, probeDist] = FindUnitLocation_mPFC(S_FNs,varargin)


% Loop through the set of provided units and identify their location in the
% brain both in a standard framework and where they sat on their respective
% probe and in absoulute DV depth into brain.
%
% Input:
% S_FNs - A cell array of the FULL file name of each unit to be localized.
%
% Outputs:
% stdLocation - An nCell x 3 Matrix of the locaiton of each unit in
%       standard space.
% netDepth - The DV distance in real space that each unit is located from
%       the surrace of brain
% peakChan - The sorted channel where each unit was localized on its
%       respective probe.
% probeDist - The distance in real space of the unit from the bottom of the
%       probe. 

% GWD March 2020

TFileRoot = '.t'; % What is the suffix on the spike time files
WaveFormSuffix = '_WF'; % What is the suffix on the WF matlab files
ChannelMapSuffix = {'ProbeA' 'ProbeC'}; % What are the Channel Map names
Probes = {'Si01' 'Si02'}; % What are the probes called in spike .t FNs; 'Si01' 'Si02'


[RootDir, DataDefName, BehavDir, TaskFileRoot, PathFileRoot, KeysFileRoot, SpikesDir, LFPDir, ChannelMapDir] = ID_GWD_Dirs;

process_varargin(varargin);

%%
dataDef = load([RootDir,DataDefName]);


% Get the Rat/SSN/WF file name for each cell
RatID = cellfun(@(x) x(1:4),S_FNs,'UniformOutput',0);
SSN = cellfun(@(x) x(1:15),S_FNs,'UniformOutput',0);

WF_FNs = cellfun(@(x,y,z) fullfile(RootDir,SpikesDir,x,y,[z TFileRoot(2:end) WaveFormSuffix]),RatID,SSN,S_FNs,'UniformOutput',0);


% Find the max depth for each rat
AllRats = unique(dataDef.RatID);
maxDepth = cell(length(AllRats),1);
for iR = 1:length(AllRats)
    dataDefIdx = ismember(dataDef.RatID,AllRats{iR});
    fullDepths = dataDef.ProbeEndDepth(dataDefIdx);
    fullDepths = cat(1,fullDepths{:});
    
    maxDepth{iR} = max(fullDepths);
end


% Initalize outputs
nCells = length(S_FNs);

stdLocation = nan(nCells,3);
netDepth = nan(nCells,1);
peakChan = nan(nCells,1);
probeDist = nan(nCells,1);

% Loop through each session and get depths for the cells in that session
uniqueSSNs = unique(SSN);
for iS = 1:length(uniqueSSNs)
    RatID = uniqueSSNs{iS}(1:4);
    cellIdx = ismember(SSN,uniqueSSNs{iS});
    dataDefIdx = ismember(dataDef.SSN,uniqueSSNs{iS});
    RatIdx = ismember(AllRats,RatID);
    
    switch dataDef.Sex{dataDefIdx}
        case 'M'
            skullThickness = .75;
        case 'F'
            skullThickness = .55;
    end
    
    chanMapFile = cell(1,length(dataDef.ProbeTarget{dataDefIdx}));
    for iP = 1:length(dataDef.ProbeTarget{dataDefIdx})
        chanMapFile{iP} = fullfile(RootDir,ChannelMapDir, [RatID '_channel_map_Si_', ChannelMapSuffix{iP},'.txt']);
    end
    
    [stdLocation(cellIdx,:), netDepth(cellIdx), peakChan(cellIdx), probeDist(cellIdx)] = localizeUnitsAlongProbe(...
        WF_FNs(cellIdx),chanMapFile,dataDef.ProbeCoordinates{dataDefIdx},maxDepth{RatIdx},...
        dataDef.ProbeDepth{dataDefIdx},'skullThickness',skullThickness,'Probes',Probes);
end


