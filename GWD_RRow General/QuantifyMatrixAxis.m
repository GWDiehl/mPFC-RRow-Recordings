
function [selectionMean, selectionStd, selectionCount, rawSelection] = QuantifyMatrixAxis(inputMtx,axis,varargin)

% Takes an input matrix and quantifies the average values along the x/y
% diagonal or perpendicular to it. Note that any non-square matrix will
% ignore the overhang on the long axis.

% GWD Feb 2021


[dimA, dimB] = size(inputMtx);

if ~isequal(dimA,dimB)
    fprintf(['Note: You did not input a square matrix so the long edge will not be fully utalized. \n',...
        ' Axis will reference with respect to the diagonal @ 1,1 \n'])
end

sqSize = min([dimA,dimB]);

chunkSize = round(sqSize/10);

centerMethod = @nanmean;
spreadMethod = @nanstd;

process_varargin(varargin);

assert(rem(chunkSize,1)==0,'Your chunk size must be an integer')

%%

sqMtx = inputMtx(1:sqSize,1:sqSize);

nBins = sqSize*2-1;

selectionMean = nan(nBins,1);
selectionStd = nan(nBins,1);
selectionCount = nan(nBins,1);

switch axis
    case {'Diag' 'Diagonal'}
        tempT = true(sqSize,chunkSize*2+1);
        selection = full(spdiags(tempT,-chunkSize:chunkSize,sqSize,sqSize));
        sqMtx(~selection) = nan;
        % Transpose so you can then take the diags
        prepData = rot90(sqMtx,1);
        
        rawSelection = nan(nBins,(chunkSize*2)^2+1);
    case {'Perp' 'Perpendicular'}
        prepData = sqMtx;        
        rawSelection = nan(nBins,sqSize*(chunkSize*2+1));
end

for iB = 1:nBins
    currDiag = iB-sqSize;
    startBin = currDiag - chunkSize;
    endBin = currDiag + chunkSize;
    selectedBins = triu(true(sqSize),startBin) & tril(true(sqSize),endBin);
    selectedData = prepData(selectedBins);
    
    validBins = ~isnan(selectedData);
    rawSelection(iB,1:sum(validBins)) = selectedData(validBins);
    
    selectionMean(iB) = centerMethod(selectedData);
    selectionStd(iB) = spreadMethod(selectedData);
    selectionCount(iB) = sum(validBins);
end

