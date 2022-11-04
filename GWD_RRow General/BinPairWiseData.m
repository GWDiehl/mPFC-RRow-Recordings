function [binnedData, pairsCount, binEdges] = BinPairWiseData(inputData,binningValues,varargin)

% Take a set of pairwise values and bin them in 2D according to a set of
% paired binning values. This was designed to look at pairwise
% relationships between cells of differing anatomical depths though it
% could be utalized more broadly.
%
% Input
% inputData - An N x 1 set of pairwise metrics you wish to bin into a
%       histogram.
% binningValues - An N x 2 set of matching values that you would like to
%       bin data acording to. Each column should correspond to values from
%       one dimention to bin against (one pair of the cell in the original
%       formulation). **Pretty sure that this should generalize out to ND
%       just fine.**
%
% Optional
% groupingMethod - annonamous function of what you would like to use to
%       group data binned together (@nanmean, @nanmedian, @nansum, etc.)
% Bins - Cell array of the bin edges you want to bin according to (passed
%       into histcn). Default uses 32 bins.
%

groupingMethod = @nanmean;
Bins = [];

process_varargin(varargin);

if isempty(Bins)
    if isempty(inputData)
        binnedData = nan(repmat(32,size(binningValues,2)));
        pairsCount = zeros(repmat(32,size(binningValues,2)));
        binEdges = binnedData;
    else
        [binnedData, binEdges] = histcn(binningValues,'AccumData', inputData, 'FUN', groupingMethod);
        pairsCount = histcn(binningValues);
    end
else
    if isempty(inputData)
        binnedData = nan(cellfun(@(x) length(x)-1,Bins));
        pairsCount = zeros(cellfun(@(x) length(x)-1,Bins));
        binEdges = binnedData;
    else        
        [binnedData, binEdges] = histcn(binningValues, Bins{:}, 'AccumData', inputData, 'FUN', groupingMethod);
        pairsCount = histcn(binningValues, Bins{:});
    end
end
