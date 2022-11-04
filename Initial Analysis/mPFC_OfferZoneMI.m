function mPFC_OfferZoneMI


%% Compute the mutual information between cell spiking and behaviors (Offer Zone)

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;

limitEventBounds = [-1 0 1 2];
nTimeSteps = [10 20 10];

overrideTSD = 'Accept';
overrideValue = [.5 1.5];
overrideDuration = [nan nan nan 1];
overrideRef = [nan nan nan 3];

nShuffles = 30; % How many shuffles of phys data do we want
shuffMethod = 'TimeShift';
commonTShift = 0;


%%% InfoTheory Specific
%Analysis Method
Method = 'PairedMI';
bins = 5;
outputDir = 'Revision\eLife\Prior\InformationTheory\';


%% Initial accept choice
currTC = 'Accept';
RRowEvent = {'OZEntry'};
EventLimits = {''};

[Results, Params] = RRow_InfoTheory_TimeWarpped('behavTC',currTC,'nBins',bins,'RRowEvent',RRowEvent,'EventLimits',EventLimits,...
    'Method',Method,'nShuffles',nShuffles,'shuffMethod',shuffMethod,'limitEventBounds',limitEventBounds,'nTimeSteps',nTimeSteps,...
    'overrideTSD',overrideTSD,'overrideValue',overrideValue,'overrideDuration',overrideDuration,'overrideRef',overrideRef);

outputFN = 'Choice_Accept';

save([RootDir,AnalysisDir,outputDir,outputFN],'Results','Params')
fprintf('\n')

%% Offer split by the choice that was made

allTCs = {'Delay' 'Value'};
RRowEvent = {'OZEntry' 'OZEntry'};
EventLimits = {'Accept' 'Skip'};

for iT = 1:length(allTCs)
    
    [Results, Params] = RRow_InfoTheory_TimeWarpped('behavTC',allTCs{iT},'nBins',bins,'RRowEvent',RRowEvent,'EventLimits',EventLimits,...
        'Method',Method,'nShuffles',nShuffles,'shuffMethod',shuffMethod,'limitEventBounds',limitEventBounds,'nTimeSteps',nTimeSteps,...
        'overrideTSD',overrideTSD,'overrideValue',overrideValue,'overrideDuration',overrideDuration,'overrideRef',overrideRef);
    
    outputFN = ['OfferByChoice_',allTCs{iT}];
    
    save([RootDir,AnalysisDir,outputDir,outputFN],'Results','Params')
    
end


%% Site Number that the rat is at

%( Deterministiclly linked to site rank, thus identical MI results)

currTC = 'SiteNumber';
RRowEvent = {'OZEntry' 'OZEntry' 'OZEntry'};
EventLimits = {'' 'Accept' 'Skip'};

[Results, Params] = RRow_InfoTheory_TimeWarpped('behavTC',currTC,'nBins',bins,'RRowEvent',RRowEvent,'EventLimits',EventLimits,...
    'Method',Method,'nShuffles',nShuffles,'shuffMethod',shuffMethod,'limitEventBounds',limitEventBounds,'nTimeSteps',nTimeSteps,...
    'overrideTSD',overrideTSD,'overrideValue',overrideValue,'overrideDuration',overrideDuration,'overrideRef',overrideRef);

outputFN = 'SiteNumber_ByChoice';

save([RootDir,AnalysisDir,outputDir,outputFN],'Results','Params')

fprintf('\n Done with Computing Choice MI \n')



