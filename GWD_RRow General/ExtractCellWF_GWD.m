
function cellWF = ExtractCellWF_GWD(cellID,varargin);


% Extract out the WF information for selected cells from the GWD mPFC data
% set. Also flips the WF of any cells that had an inverted WF so that the
% peak/valley/rise and fall times are straightforward.

% GWD Nov 2021

[RootDir, ~, ~, ~, ~, ~, SpikesDir] = ID_GWD_Dirs;
TFileRoot = '.t';
WaveFormSuffix = '_WF';

WFPeak = [4 12];

process_varargin(varargin);

%%

nCells = length(cellID);
cellWF = cell(nCells,1);

% Get the Rat/SSN/WF file name for each cell
RatID = cellfun(@(x) x(1:4),cellID,'UniformOutput',0);
SSN = cellfun(@(x) x(1:15),cellID,'UniformOutput',0);
WF_FNs = cellfun(@(x,y,z) fullfile(RootDir,SpikesDir,x,y,...
    [z TFileRoot(2:end) WaveFormSuffix]),RatID,SSN,cellID,'UniformOutput',0);

for iC = 1:nCells
    fullWF = load(WF_FNs{iC});
    [~, peakChan] = max(max(abs(fullWF.Mean)));
    tempWF = fullWF.Mean(:,peakChan);
    
    % Negative dir WF
    negPeak = min(tempWF(WFPeak(1):WFPeak(2)));
    % Positive dir WF
    posPeak = max(tempWF(WFPeak(1):WFPeak(2)));
    
    % PosDir, flip sign
    if posPeak > -negPeak
        tempWF = -tempWF;
    end
    
    % Double check that the peak is in the right place
    [~, peakIdx] = min(tempWF);
    if peakIdx < WFPeak(1) || peakIdx > WFPeak(2)
        tempWF(:) = nan;
    end
    cellWF{iC} = tempWF;
end

