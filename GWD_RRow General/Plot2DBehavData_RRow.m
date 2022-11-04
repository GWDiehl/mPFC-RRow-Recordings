

function [fh,fullData,groupingID,grpResults,binEdges] = Plot2DBehavData_RRow(variables,varargin)


% GWD March 2021

assert(iscell(variables) & length(variables) <4,...
    'You need 1-3 variables, entered as strings in a cell')

 % Positive Next lap; Negative Previous lap
taskStep = zeros(size(variables));

groupingFunc = @nanmean;

binEdges = cell(1,length(variables));

groupSplit = 'ByRat';
plotGrps = 1;
normalizeResults = 'None';


%% Select and load the requisite behavioral data set

% Prelim Directoires & DataDef
[rootDir, DataDefName,BehavDir] = ID_GWD_Dirs;
dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

SSNs = dataDef.SSN(OI);

% Get any behavioral Task Data that may be needed
[FullLapData, FullSessionData, FullIdPhi] = collectRRowDataToAnalyze('restrictionField',{'SSN'},'restrictionValue',{SSNs});

process_varargin(varargin);


%%


fullData = cell(1,length(variables));


switch groupSplit
    case 'ByRat'
        allRats = unique(FullSessionData.RatID,'stable');
        groupingID = cellfun(@(x) find(ismember(allRats,x)),FullSessionData.RatID);
    case 'BySess'
        groupingID = [1:length(FullSessionData.RatID)]';
    case 'None'
        groupingID = ones(length(FullSessionData.RatID),1);
    otherwise
        error('Split case not implemented')
end

nGroups = max(groupingID);
plotGrpColors = jet(nGroups);

if nGroups > 15
    plotGrps = 0;
end

for iV = 1:length(variables)
    
    [currVariable,defaultEdges] = ExtractBehavDataFromFull(FullLapData,FullSessionData,variables{iV});
    
    sessTime = ExtractBehavDataFromFull(FullLapData,[],'SessionTime');
%     currVariable(sessTime<600) = nan;
    currVariable(sessTime>600) = nan;
    
    if taskStep(iV)
        if taskStep(iV) > 0
            prefix = 'Next';
        elseif taskStep(iV) < 0
            prefix = 'Prev';
        end
        variables{iV} = sprintf('%s_%s%d',variables{iV},prefix,abs(taskStep(iV)));
        currVariable = advanceRRowData(currVariable,-taskStep(iV));
    end
    
    fullData{iV} = currVariable;
    if isempty(binEdges{iV})
        binEdges{iV} = defaultEdges;
    end    
end


nBins = cellfun(@(x) length(computeBinCenters(x)),binEdges);
grpResults = nan([nGroups nBins(2:end)]);

for iG = 1:nGroups
    currIdx = groupingID == iG;
    
    dependData = reshape(fullData{1}(currIdx,:),[],1);
    
    independData = [];
    for iV = 2:length(variables)
        independData = cat(2,independData,reshape(fullData{iV}(currIdx,:),[],1));
    end
    
    validIdx = ~isnan(sum(cat(2,independData,dependData),2));
    
    temp = histcn(independData(validIdx,:), binEdges{2:end}, 'AccumData', dependData(validIdx), 'FUN',groupingFunc);
    counts = histcn(independData(validIdx,:), binEdges{2:end});
    temp(counts==0) = nan;
    switch length(variables)
        case 2
            temp = temp(1:nBins(2));
            grpResults(iG,:) = temp;
        case 3
            temp = temp(1:nBins(2),1:nBins(3));
    grpResults(iG,:,:) = temp;
    end
    
    switch normalizeResults
        case 'SessCount'
            grpResults(iG,:) = grpResults(iG,:)/sum(currIdx);
        case 'LapCount'
            grpResults(iG,:) = grpResults(iG,:)/sum(validIdx);
    end
end

dataMean = nanmean(grpResults,1);
dataSEM = nanstd(grpResults,[],1)./sqrt(sum(~isnan(grpResults),1));

fh = figure;
if length(variables) == 2
    if plotGrps
        hold on
        for iG = 1:nGroups
            plot(computeBinCenters(binEdges{2}),grpResults(iG,:),'color',plotGrpColors(iG,:))
        end
    end
    errorbar(computeBinCenters(binEdges{2}),dataMean,dataSEM,'k','LineWidth',1.5)
    xlabel(variables{2})
    ylabel(variables{1})
else
    s = imagesc(computeBinCenters(binEdges{2}),computeBinCenters(binEdges{3}),squeeze(dataMean)');
    set(s,'AlphaData',~isnan(squeeze(dataMean)'));
    
    title(variables{1})
    xlabel(variables{2})
    ylabel(variables{3})
    axis xy
end
