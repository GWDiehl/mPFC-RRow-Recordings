function [StandardLocation, netDepth, peakChan, probeDist] = localizeUnitsAlongProbe(WF_FNs,channelMap,ProbesCoordinates,EndDepths,CurrentDepths,varargin)

% Find the anatomical location of each recorded unit based on the largest
% WF across a set of probe channels. Take the absolute depth of this
% recording and convert it into a location in a standariizded 3D space.
%
% INPUT
% WF_FNs - Cell array of the file names of the WF calcualtions. Expects a
%       structure with .Mean field of Time x Channels for each unit.
%
% ** Each of the following are expected to be cell arrays where each
% element corresponds to a given probe **
% channelMap - File name of the Intan style channel map.
% ProbesCoordinates - 3 x 2 matrix of the 3D coordinates of the top(:,1)
%       and bottom(:,2) of the probe in standard space
% **These two would be vectors of length(nProbes)**
% EndDepths - Final depth of the probe in real space (in mm)
% CurrentDepths - Current depth of the bottom of the probe in real space
%       (in mm)
%
% OPTIONAL INPUT
% skullThickness - Our measurements tend to be from the top of skull
%       whereas the top of the probe reflects the top of brain. Subtract
%       skull thickness from real space distances
% Probes - Cell array of the FileName parts corresponding to each probe to
%       be fond in the file names
%
% OUTPUTS
% StandardLocation - An nCell x 3 Matrix of the locaiton of each unit in
%       standard space.
% netDepth - The DV distance in real space that each unit is located from
%       the surrace of brain
% peakChan - The sorted channel where each unit was localized on its
%       respective probe.
% probeDist - The distance in real space of the unit from the bottom of the
%       probe.

% GWDiehl March 2020


skullThickness = .75;
Probes = {'Si01' 'Si02'}; %'Si01' 'Si02'

process_varargin(varargin);

%%

nCells = length(WF_FNs);

% If things came from multiple probes but only a single channel map,
% use empty to group them all from the same thing.
if isempty(Probes)
    probeIdx = ones(nCells,1);
else
    % Identifiy which sensor and target each cell corresponds to
    probeIdx = nan(nCells,1);
    for iP = 1:length(Probes)
        probeIdx(contains(WF_FNs,Probes{iP})) = iP;
    end
end

% Grab sorted channel maps, and distances for each site on each probe
usedProbes = length(channelMap);
SortChanNum = cell(1,usedProbes); % Sorting of channel numbers
ChanFromTop = cell(1,usedProbes); % Depth of each sorted channel from top of probe
ChanFromBottom = cell(1,usedProbes); % From the bottom of the probe going up
for iP = 1:usedProbes    
    [~, channel_map, ~, ~, ~, ~, ycoords] = GetChannelMap_kilosort(channelMap{iP});
    
    [sortedDepth, chanOrder] = sort(ycoords);
    SortChanNum{iP} = channel_map(chanOrder)+1-min(channel_map);
    ChanFromTop{iP} = sortedDepth/1000; % Convert to mm
    ChanFromBottom{iP} = ChanFromTop{iP}(end) - ChanFromTop{iP};
end


% Localize each unit to its respective site based on WF location
peakChan = nan(nCells,1); % Sorted channel with largest WF
probeDist = nan(nCells,1); % Distance from the bottom of the probe
netDepth = nan(nCells,1); % Net depth of the unit from the brain surface (skull subtracted)
StandardLocation = nan(nCells,3); % 3D coordinate in standardized space based on scaling along the trajectory
for iC = 1:nCells
    WF = load(WF_FNs{iC});
    sortedWF = WF.Mean(:,SortChanNum{probeIdx(iC)});
    [~, peakChan(iC)] = max(max(abs(sortedWF),[],1));
    
    % Distance from the bottom of the probe and current depth of the unit;
    % Unit measures
    probeDist(iC) = ChanFromBottom{probeIdx(iC)}(peakChan(iC));
    netDepth(iC) = CurrentDepths(probeIdx(iC)) - probeDist(iC) - skullThickness;
    
    
    % Top, bottom, final depth of the probe; Probe measures
    refProbeTop = ProbesCoordinates{probeIdx(iC)}(:,1);
    refProbeBottom = ProbesCoordinates{probeIdx(iC)}(:,2);
    refProbeDepth = EndDepths(probeIdx(iC)) - skullThickness; % Real Space
    
    % Find the location in standardized 3D space
    StandardLocation(iC,:) = rescaleIn3D(netDepth(iC),refProbeDepth,refProbeTop,refProbeBottom);
end
