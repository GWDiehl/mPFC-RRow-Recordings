function [avgThresh, threshSpread, threshCorr] = PlotThresholdMetrics(Thresholds,groupIDs)

% Plots three different sets of analyses regarding RRow Thresholds. Plots
% the distribution + Mean/Std for across sessions for each rat (unique
% groupIDs), The varience in threholds across day or RR, The correlation
% between Threshold sets from day to day, split between within the same rat
% (expected similar profiles) and across different rats (expect ~chance
% level).


correlationType = 'Pearson';

correlationData = corr(Thresholds','type',correlationType);
correlationData(diag(true(1,size(Thresholds,1)))) = nan;

nSites = size(Thresholds,2);
nGroups = max(groupIDs);
grpColors = jet(nGroups);

avgThresh = figure; hold on
threshSpread = figure;hold on
threshCorr = figure;hold on

stepSize = 1/nGroups/2;
offSet = stepSize*nGroups/2;

sameRat = false(length(groupIDs));
xRat = true(length(groupIDs));

fullSessVar = nan(1,nGroups);
fullFlavorVar = nan(1,nGroups);

for iG = 1:nGroups
    currIdx = groupIDs == iG;
    currData = Thresholds(currIdx,:);
    
    figure(avgThresh)     
    plot([1:nSites]-offSet+stepSize*(iG-.5),currData,'.','color',grpColors(iG,:))
    errorbar([1:nSites]-offSet+stepSize*(iG-.5),nanmean(currData,1),nanstd(currData,[],1)/sqrt(size(currData,1)),'color',grpColors(iG,:))
    
    figure(threshSpread)
    sessVar = nanvar(currData,[],1);
    sessError = nanstd(sessVar)/sqrt(length(sessVar));
    flavorVar = nanvar(currData,[],2);
    flavorError = nanstd(flavorVar)/sqrt(length(flavorVar));
    errorbar(nanmean(sessVar),nanmean(flavorVar),flavorError,flavorError,sessError,sessError,'color',grpColors(iG,:))
    fullSessVar(iG) = nanmean(sessVar);
    fullFlavorVar(iG) = nanmean(flavorVar);
    
    sameRat(currIdx,currIdx) = 1;
    xRat(currIdx,currIdx) = 0;
end

[p, h, stats] = signrank(fullSessVar,fullFlavorVar);

figure(avgThresh)
xlabel('Restaurant Num')
ylabel('Threshold')

figure(threshSpread)
xlabel('Var across Sess')
ylabel('Var across RR')
plotIDLine({'k'})


figure(threshCorr)
histogram(correlationData(sameRat),-1:.1:1,'Normalization','probability','DisplayStyle','stairs')
histogram(correlationData(xRat),-1:.1:1,'Normalization','probability','DisplayStyle','stairs')
title('Day-To-Day Correlation of Threshold profiles')
legend({'Same Rat', 'Cross Rat'},'location','best')
xlabel('Correlation between Threshold profiles')
ylabel('Proportion')
