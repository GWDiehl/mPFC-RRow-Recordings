function binCenters = computeBinCenters(binEdges)

% Take bin edges (vector) and compute the bin centers

binCenters = binEdges(2:end) - diff(binEdges)/2;