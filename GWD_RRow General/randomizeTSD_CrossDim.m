function outputTSD = randomizeTSD_CrossDim(inputTSD)

% Takes an input TSD and does a full randomize permuation of ALL data
% entries across ALL dimensions. The same data are used, such that the same
% overall probability distributions are maintained but matching data
% entries in the original will no longer be matching in the output.
%
% All new entries are taken from the original data and nan entries remain
% unchanged. Data are randomized across all dimensions.
% 
% See also 'shuffleTSD' for a shuffling of data values while maintining
% matching enteres.

% GWD March 2022

inputData = inputTSD.data;

validEntries = ~isnan(inputData(:));

origData = inputData(validEntries);
inputData(validEntries) = randsample(origData,length(origData));

outputTSD = tsd(inputTSD.range,inputData);