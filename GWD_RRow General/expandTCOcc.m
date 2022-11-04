function expandedOcc = expandTCOcc(basicOcc,nCells)

% Take a single iteration of occupacy output from TuningCurve functions and
% replicate it according to number of cells to match the output
% SpikeCounts/AvgRate. Can also pass in a cell array of data corresponding
% to multiple series (sessions, TCs, etc.)

% GWD Oct 2020


if iscell(basicOcc(1))
    nSess = length(basicOcc);
    expandedOcc = cell(nSess,1);
    for iS = 1:nSess
        temp = basicOcc{iS};
        expandedOcc{iS} = reshape(repmat(temp(:),1,nCells(iS))',[nCells(iS),size(temp)]);
    end
else
    expandedOcc = reshape(repmat(basicOcc(:),1,nCells)',[nCells,size(basicOcc)]);
end