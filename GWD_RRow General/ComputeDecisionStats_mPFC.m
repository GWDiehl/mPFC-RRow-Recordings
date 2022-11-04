function [rmANOVA,pVals,sigTime,ANOVA,ANOVA_stats,multiCompare] = ComputeDecisionStats_mPFC(mainData,groupingFields,selectedWindow,varargin)

alphaVal = 0.05;
multiCorrect = 1;
tailedTest = 'right';

process_varargin(varargin);


% Collect data for stats
nCells = cellfun(@(x) size(x,1),mainData);
nGroups = length(groupingFields);
nTimes = size(mainData{1},2);

pVals = nan(nGroups,nTimes);
significant = false(nGroups,nTimes);
if multiCorrect
    alphaVal = alphaVal/nTimes; % MultiCompare adjusted significnace Val
end

fullDataList = arrayfun(@(x) repmat(groupingFields(x),1,nCells(x))',1:nGroups,'UniformOutput',0);
fullDataList = cat(1,fullDataList{:});

fullOutputData = cat(1,mainData{:});


% Repeated Measures ANOVA
dataTable = cell2table(cat(2,fullDataList,num2cell(fullOutputData)),...
    'VariableNames',[{'Subgroups'},arrayfun(@(x) sprintf('Time%d',x),[1:nTimes],'UniformOutput',0)]);

rm = fitrm(dataTable,sprintf('Time1-Time%d~Subgroups',nTimes),'WithinDesign',[1:nTimes]');
rmANOVA = ranova(rm);

for iG = 1:nGroups
    for iT = 1:nTimes
        pVals(iG,iT) = signrank(mainData{iG}(:,iT),0,'tail',tailedTest);
    end
end
sigTime = pVals < alphaVal;

% Standard 1-Way ANOVA computed over the decision period Area Under Curve
selectedData = nansum(fullOutputData(:,selectedWindow),2);
[p,ANOVA,ANOVA_stats] = anova1(selectedData,fullDataList,'off');
multiCompare = multcompare(ANOVA_stats,'display','off');

[~,~,idx] = unique(fullDataList);
[p,~,stats] = ranksum(selectedData(idx==1),selectedData(idx==2));