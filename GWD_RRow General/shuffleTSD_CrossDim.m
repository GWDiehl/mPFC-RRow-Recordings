function outputTSD = shuffleTSD_CrossDim(inputTSD)

% Takes an input TSD and shuffles identity of the data across ALL
% dimensions using a consistant shuffle pattern. All original examples of X
% will now be Y. Each dimension is treated independently (i.e. each cell is
% shuffled internally).
%
% All new entries are taken from the original data and nan entries remain
% unchanged. Data are taken across ALL dimensions
%
% See also 'randomizeTSD' for a full random permuataion of all data
% entries.

% GWD March 2022

outputTSD = inputTSD;
outputTSD.D(:) = nan;

inputData = inputTSD.data;

validEntries = ~isnan(inputData(:));

origVal = unique(inputData(validEntries,iV));
randVal = randsample(origVal,length(origVal));
for iE = 1:length(randVal)
    outputTSD.D(inputData(:) == origVal(iE)) = randVal(iE);
end