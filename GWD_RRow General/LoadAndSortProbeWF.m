
function [sortedWFs,peakChan] = LoadAndSortProbeWF(WF_FNs,chanMap_FNs,varargin)

% For an array of spike WF file names and corresponding channel maps, load
% the WF data and sort according to the channel map listing. Also extract
% out the peak channel of the sorted data

% GWD March 2021


expandOutput = 0;
process_varargin(varargin);

%%

nCells = length(WF_FNs);

sortedWFs = cell(nCells,1);
peakChan = nan(nCells,1);

for iC = 1:nCells    
    [~, channel_map, ~, ~, ~, ~, ycoords] = GetChannelMap_kilosort(chanMap_FNs{iC});
    
    % Sort by depth
    [~, chanOrder] = sort(ycoords);
    % Sort the channel numbers and convert them to Idx values
    sortChanNum = channel_map(chanOrder)+1-min(channel_map);
    
    WFData = load(WF_FNs{iC});
    sortedWFs{iC} = WFData.Mean(:,sortChanNum);
    [~, peakChan(iC)] = max(max(abs(sortedWFs{iC}),[],1));
end

if nCells == 1 && expandOutput
   sortedWFs = sortedWFs{1};
end