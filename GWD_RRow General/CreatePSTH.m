function binCounts = CreatePSTH(spikeData,stimTimes,histEdges)


nCells = length(spikeData);
nEvents = length(stimTimes);
nBins = length(histEdges)-1;
binCounts = zeros(1,nBins);

for e = 1:nEvents
    spikeData = cellfun(@(x) restrict(x,stimTimes(e)+histEdges(1),stimTimes(e)+histEdges(end)),spikeData);
    
    for c = 1:nCells
        eventBins = histcounts(spikeData(c).T,histEdges);
        binCounts = binCounts + eventBins;
    end
end