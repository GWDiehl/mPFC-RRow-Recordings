function outputTSD = randomizeTSD(inputTSD)

% Takes an input TSD and does a randomize permuation of all data entries
% for each dimension (i.e. each cell is randomized independently). The same
% data are used, such that the same overall probability distributions are
% maintained but matching data entries in the original will no longer be
% matching in the output.
%
% All new entries are taken from the original data and nan entries remain
% unchanged. Each TSD dimension is handled independently.
% 
% See also 'shuffleTSD' for a shuffling of data values while maintining
% matching enteres.

% GWD July 2020

inputData = inputTSD.data;
nVar = size(inputData,2);

for iV = 1:nVar    
    validEntries = ~isnan(inputData(:,iV));
    
    origData = inputData(validEntries,iV);
    inputData(validEntries,iV) = randsample(origData,length(origData));
end

outputTSD = tsd(inputTSD.range,inputData);