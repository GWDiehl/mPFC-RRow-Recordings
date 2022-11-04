

function [FRCorrDT,minPairwiseOcc,FieldID,ResultID,OccID,VarID] = ResultsFieldNames(dataType)

% Spit out the standard field name info for Results structures for various
% data types.

% GWD Feb 2021


FRCorrDT = nan;
minPairwiseOcc = nan;
FieldID = '';
ResultID = '';
OccID = '';
VarID = '';

switch dataType
    case {'FR' 'FRs'}
        FRCorrDT = .05;
    case {'TC' 'TCs' 'TuningCurves'}
        minPairwiseOcc = 1;
        FieldID = 'FieldNames';
        ResultID = 'AvgRate';
        OccID = 'Occ';
        VarID = 'RateVar';
    case {'PETH' 'PETHs' 'PETH_TW'}
        minPairwiseOcc = .3;
        FieldID = 'FullEventID';
        ResultID = 'AvgRate';
        OccID = 'AvgOcc';
        VarID = 'RateVar';
    case{'MI' 'MutualInfo' 'MI_TW'}
        minPairwiseOcc = .3;
        FieldID = 'FullEventID';
        ResultID = 'Info';
        OccID = 'AvgOcc';
    case {'Corr' 'Correlation'}
        minPairwiseOcc = .15;
        FieldID = 'FullEventID';
        ResultID = 'AvgCorr';
        OccID = 'AvgOcc';
    case 'General'
        minPairwiseOcc = .3;
        FieldID = 'FullEventID';
        ResultID = 'AvgResults';
        OccID = 'AvgOcc';    
        VarID = 'ResultsVar';
    otherwise
        error('Unknown Data Type')
end

if contains(dataType,'TW')
    OccID = '';
end