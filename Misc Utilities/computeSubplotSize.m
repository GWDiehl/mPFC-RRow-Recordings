function [subplotA, subplotB] = computeSubplotSize(nElements)

% Takes a number of plotted figures and computes a logical breakdown of
% subplotting

% GWD Dec 2020


breakpoint = ceil(sqrt(nElements));
subplotA = ceil(nElements/breakpoint);
subplotB = ceil(nElements/subplotA);