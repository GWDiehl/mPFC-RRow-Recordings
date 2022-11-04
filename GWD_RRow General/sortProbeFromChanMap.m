function [sortedChanNum, sensorNum, depth] = sortProbeFromChanMap(chanMap,varargin)

ChannelsToSort = [];
respectSensors = true;
process_varargin(varargin);


[~, channel_map, sensor, ~, ~, ~, ycoords] = GetChannelMap_kilosort(chanMap);
nChan = length(channel_map);
if ~isempty(ChannelsToSort)
    channel_map = channel_map(ChannelsToSort);
    sensor = sensor(ChannelsToSort);
    ycoords = ycoords(ChannelsToSort);
end

sortedChanNum = nan(nChan,1);
sensorNum = nan(nChan,1);
depth = nan(nChan,1);

if respectSensors
    allSensors = unique(sensor);
    for iS = 1:length(allSensors)
        toSort = sensor == allSensors(iS);
        tempY = ycoords(toSort);
        tempChanMap = channel_map(toSort);
        
        [sortedDepth, chanOrder] = sort(tempY);
        chanNum = tempChanMap(chanOrder)+1-min(channel_map);
        
        sortedChanNum(toSort) = chanNum;
        sensorNum(toSort) = allSensors(iS);
        depth(toSort) = sortedDepth;
    end
else    
    [sortedDepth, chanOrder] = sort(ycoords);
    chanNum = channel_map(chanOrder)+1-min(channel_map);
    
    sortedChanNum = chanNum;
    sensorNum = sensor(chanOrder);
    depth = sortedDepth;
end
