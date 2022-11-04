function [selectedData,fieldName,outputScale] = ExtractRegressionOutputs(results,comparison,metric)

switch metric
    case 'InModel'
        switch comparison
            case 'Stepwise'
                selectedData = cat(1,results.(comparison).InModel{:});
            case {'OrderedStep' 'Individual'}
                selectedData = cat(1,results.(comparison).pVal{:})<0.05;
        end
        fieldName = 'Prop In Model';
        outputScale = [0 1];
    case {'RankOrder' 'RankOrder_InModel'}
        modelPVal = cat(1,results.(comparison).pVal{:});
        selectedData = nan(size(modelPVal));
        for iX = 1:size(rankOrder,1)
            [~,order] = sort(modelPVal(iX,:));
            selectedData(iX,order) = flip(1:nBehavVar);
        end
        if strcmp(fieldOfInterest,'RankOrder_InModel')
            selectedData(modelPVal>0.05) = notInModelVal;
        end
        fieldName = 'Rank Order';
        outputScale = [0 nBehav];
    case 'Correlation'
        selectedData = cat(1,results.(comparison).Correlation{:});
        fieldName = 'Correlation';
        outputScale = [-.3 .3];
    case 'AbsCorrelation'
        selectedData = abs(cat(1,results.(comparison).Correlation{:}));
        fieldName = 'Abs Correlation';
        outputScale = [0 .5];
    case 'ExplainedVar'
        selectedData = cat(1,results.(comparison).Correlation{:}).^2;
        fieldName = 'Explained Var';
        outputScale = [0 .1];
    case 'MutualInfo'        
        selectedData = cat(1,results.(comparison).MutualInfo{:});
        fieldName = 'Mutual Info';
        outputScale = [0 .2];
end