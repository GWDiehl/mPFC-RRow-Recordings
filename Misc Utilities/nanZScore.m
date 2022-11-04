function zscoreOut = nanZScore(inputVector)


% GWD March 2021

assert(isvector(inputVector),['nanZScore requires a simple vector.'...
    'If you want to deal with a matrix do it byhand'])

validEntries = ~isnan(inputVector);
if sum(validEntries) == 0
    zscoreOut = nan;
    return
end
zscoreOut = nan(size(inputVector));

zscoreOut(validEntries) = zscore(inputVector(validEntries));
