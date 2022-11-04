
function [groupedTEVals,fh_byRegion,fh_byDistance,ANOVA,ANOVA_Stats,multiCompare] = PlotTEAcrossSubregions(pairwiseVals,cellPair,cellID,anatGroups,varargin)


xRange = [-1e-4 6e-4];
grpColors = loadGWD_RRowColormap;

% Give extra spacing in your jet color scheme as it seems to help
% descriminability
jetSpacing = 2;

process_varargin(varargin);

%%


groupingID = groupPFCCellsByDepth(cellID,anatGroups);
nAnat = length(anatGroups);
[subA,subB] = computeSubplotSize(nAnat);

groupedTEVals = cell(nAnat);
for iA = 1:nAnat
    validPair_A = ismember(cellPair(:,1),find(groupingID(:,iA)));
    for iB = 1:nAnat
        validPair_B = ismember(cellPair(:,2),find(groupingID(:,iB)));
        groupedTEVals{iA,iB} = pairwiseVals(validPair_A & validPair_B);
    end
end

fh_byRegion = figure;
for iB = 1:nAnat
    subplot(subA,subB,iB); hold on
    for iA = 1:nAnat
        [f1,x1] = ecdf(groupedTEVals{iA,iB});
        plot(x1,f1,'color',grpColors.(anatGroups{iA}))
    end
    
    title(sprintf('%s Info Transfer to:',anatGroups{iB}))
    
    xlim(xRange)
    plotVertLine(0,{'--k'})
    plotHorizLine(0.5,{'--k'})
    
    if iB == 1
        legend(anatGroups,'location','best')
    end
    xlabel('Normalized Transfer Entropy')
    ylabel('Cummulative Proportion')
end



intraRegion = diag(true(nAnat,1));
entryNum = reshape(1:nAnat^2,nAnat,nAnat);
colorSets = jet([nAnat*jetSpacing]-1);
colorSets = colorSets([1:jetSpacing:size(colorSets,1)],:);

fh_byDistance = figure;
subplot(1,2,1); hold on

fullPairData = [];
groupingID = [];

tempData = cat(1,groupedTEVals{intraRegion});
[f1,x1] = ecdf(tempData); plot(x1,f1);
fullPairData = cat(1,fullPairData,tempData);
groupingID = cat(1,groupingID,repmat(1,length(tempData),1));


tempData = cat(1,groupedTEVals{~intraRegion});
[f1,x1] = ecdf(tempData); plot(x1,f1);
fullPairData = cat(1,fullPairData,tempData);
groupingID = cat(1,groupingID,repmat(2,length(tempData),1));

[p,ANOVA,ANOVA_stats] = anova1(fullPairData,groupingID,'off');
multiCompare = multcompare(ANOVA_stats,'display','off');

[p, h, stats] = ranksum(fullPairData(groupingID==1),fullPairData(groupingID==2));


xlim(xRange)
plotVertLine(0,{'--k'})
plotHorizLine(0.5,{'--k'})
legend({'Within Region' 'Across Region'},'location','best')
xlabel('Normalized Transfer Entropy')
ylabel('Cummulative Proportion')

fullPairData = [];
groupingID = [];
subplot(1,2,2); hold on
for iD = 1:nAnat
    selectedIdx = unique(cat(1,diag(entryNum,iD-1),diag(entryNum,1-iD)));    
    tempData = cat(1,groupedTEVals{selectedIdx});
    
    [f1,x1] = ecdf(tempData);
    plot(x1,f1,'color',colorSets(iD,:))
    
    fullPairData = cat(1,fullPairData,tempData);
    groupingID = cat(1,groupingID,repmat(iD-1,length(tempData),1));
end

[p,ANOVA,ANOVA_stats] = anova1(fullPairData,groupingID,'off');
multiCompare = multcompare(ANOVA_stats,'display','off');

[r,p] = nancorr(fullPairData,groupingID);

[f1,x1] = ecdf(cat(1,groupedTEVals{~intraRegion}));
plot(x1,f1,'--k','linewidth',1)

xlim(xRange)
plotVertLine(0,{'--k'})
plotHorizLine(0.5,{'--k'})
legend(cat(2,arrayfun(@(x) sprintf('%d Subregions away',x-1),1:nAnat,'UniformOutput',0),'All Across'),'location','best')
xlabel('Normalized Transfer Entropy')
ylabel('Cummulative Proportion')





