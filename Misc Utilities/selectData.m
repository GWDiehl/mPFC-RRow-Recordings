function data = selectData(originalData,selection,varargin)

% Selects a set of data within a matrix of 'originalData' based on a
% logical matrix 'selection'. In doing so it will maintain the
% size/structure of the original data input replacing any unselected
% entries with nan.

% GWD 2019

dummyValue = nan;

process_varargin(varargin);

data = originalData;
data(~selection) = dummyValue;