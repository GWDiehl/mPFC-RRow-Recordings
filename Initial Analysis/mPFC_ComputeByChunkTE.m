function mPFC_ComputeByChunkTE

%% Transfer entropy segregated by maze chunk

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;

Method = 'PairedMI';
nShuffles = 30;
shuffMethod = 'TimeShift';

dt = .01;

outputRoot = 'Revision\eLife\TransferEntropy\MazeChunks\';

% Make sure the directory exists
if ~exist([RootDir,AnalysisDir,outputRoot],'dir')
    mkdir([RootDir,AnalysisDir,outputRoot])
end

RRowEvents = {'All' 'OfferZone' 'WaitZone' 'LingerZone' 'AdvanceZone'}; %  
EventLimits = {'' 'Skip' 'Accept' 'Earn' 'Quit'};

for iE = 1:length(RRowEvents)
    
    currEvents = repmat(RRowEvents(iE),1,length(EventLimits));
    
    [Results, Params] = test_FullZoneTE('RRowEvent',currEvents,'EventLimits',EventLimits,...
        'nShuffles',nShuffles,'shuffMethod',shuffMethod,'dt',dt,'validCells',{'No_VO'});
    
    currFN = ['PairwiseTE_FullEvent_',RRowEvents{iE}];
    
    save([RootDir,AnalysisDir,outputRoot,currFN],'Results','Params')
    fprintf('\n Done with %s \n\n',RRowEvents{iE})
end
