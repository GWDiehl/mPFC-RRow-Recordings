function outputTSD = shuffleTSD(inputTSD)

% Takes an input TSD and shuffles identity of the data for each dimension using a
% consistant shuffle pattern. All original examples of X will now be Y.
% Each dimension is treated independently (i.e. each cell is shuffled
% internally).
%
% All new entries are taken from the original data and nan entries remain
% unchanged. Each TSD dimension is handled independently.
%
% See also 'randomizeTSD' for a full random permuataion of all data
% entries.

% GWD July 2020

outputTSD = inputTSD;
outputTSD.D(:) = nan;

inputData = inputTSD.data;
nVar = size(inputData,2);

for iV = 1:nVar
    validEntries = ~isnan(inputData(:,iV));
    
    origVal = unique(inputData(validEntries,iV));
    randVal = randsample(origVal,length(origVal));
    for iE = 1:length(randVal)        
        outputTSD.D(inputData(:,iV) == origVal(iE),iV) = randVal(iE);
    end
end