function [FeederRank, RankOrder] = ExtractRRowRankOrdering(SortingData,varargin)


% Input a matrix of behavioral data to be used for rank ordering RRow
% feeder sites. This is typically the threshold, number of rewards earned,
% number of site visits/accepts, etc. Matrix should be nSess x nSites x
% nSortingVariables. Sorting can be in descending (default) or ascending
% order based on varargin, but the same direction must be used for all
% sorting variables. Sessions can also be grouped with an input variable
% 'grouping'.
%
% Outputs the ranking (1: Best, N: Worst) of each original site number and
% the listing of original sit numbers in ranked order from best to worst.

% GWD 2021

grouping = 1:size(SortingData,1);
direction = 'descend';
process_varargin(varargin);

[nEntries,nSites,nVariables] = size(SortingData);

FeederRank = nan(nEntries,nSites);
RankOrder = nan(nEntries,nSites);

groupIDs = unique(grouping);
for g = 1:length(groupIDs)
    if iscell(grouping(1))
        groupMembers = ismember(grouping,groupIDs{g});
    elseif isnumeric(grouping(1))
        groupMembers = grouping == groupIDs(g);
    end
    
    groupData = squeeze(nanmean(SortingData(groupMembers,:,:),1));
    if nVariables == 1
        [~, groupOrder] = sort(groupData',direction);
    else
        [~, groupOrder] = sortrows(groupData,direction);
    end
    groupRank = arrayfun(@(x) find(x == groupOrder),1:nSites);
    
    RankOrder(groupMembers,:) = repmat(groupOrder',sum(groupMembers),1);
    FeederRank(groupMembers,:) = repmat(groupRank,sum(groupMembers),1);
end