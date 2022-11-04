function [beta,corrVal,pVal,residuals] = OrderedStepwiseRegression(xData,yData,modelOrder,varargin)


% Perform a pseudo stepwise regression in which predictor variables are
% regressed against a set of response values in a perscribed order. For
% each predictor variable data are correlated (linear) to yield a
% correlation and pVal. They are then regressed according to the polynomial
% order indicated to produce the regression beta weights. The prediction is
% then subtracted from the used responses to elicit the remaining
% residuals. These residuals are then used in the subsequent correlations
% and regressions.
%
% INPUTS
% xData - an nSamples x mVariables matrix of the predictor (independent
%       variable) data
% yData - an nSamples vector of the response (dependent variable) data
% modelOrder - The polynomial order of the regression to be used. Note that
%       this will apply to the regression and the subtraction to produce
%       residuals, but the correlation and accompanying pVal will also be a
%       simple linear correlation
%
% OPTIONAL INPUTS
% entriesToUse - A logical array of the entries to be USED in the fitting
%       of the correlation/regression. Note that all of the entred yData is
%       regressed out but some of these could be omitted in the fitting of
%       the regression.
%
%
% OUTPUTS 
% beta - The beta weights of the fit regression model for each
%       predictor varaible. Note this will be order+1 entires in decreasing
%       order as output by matlab polyfit
% corrVal - The correlation for each predictor according to the remaining
%       residuals in the response data
% pVal - The pVals for each corresponding correlation
% residuals - The residual values at each stage of the stepwise regression
%       process. Note that these are the residual AFTER regressing out the
%       particular variable of interest.

% GWD Cinco de Mayo! 2022

% Maximum beta weight that can be fit in the polynomial regression. This
% should be a BIG number that will effectivly catch ill fit relations. If
% you REALLY REALLY want to hold BIG fits then set this to inf as an input
% arg.
maxPolynomialBeta = 1e14;

entriesToUse = [];
useZTransform = 0;

process_varargin(varargin);


%%

nVar = size(xData,2);
nSamples = length(yData);

assert(isequal(size(xData,1),length(yData)),'Your predictor and response data must be the same length')

beta = nan(nVar,modelOrder+1);
corrVal = nan(nVar,1);
pVal = nan(nVar,1);
residuals = nan(nSamples,nVar);

if isempty(entriesToUse)
    entriesToUse = true(nSamples,1);
else
    assert(length(entriesToUse) == nSamples,'Your entries to use needs to be a logical that matches the number of samples in your data')
end

% zTransform the xData so that things will be well conditioned/not run into
% issues with machine precision
if useZTransform
    meanX = nanmean(xData,2);
    stdX = nanstd(xData,[],2);
    xData = (xData - meanX)./stdX;
end

for iR = 1:nVar
    origYData = yData;
    validSamples = ~isnan(xData(:,iR)+yData);
    
    yData(~validSamples) = nan;
    [corrVal(iR),pVal(iR)] = corr(xData(validSamples&entriesToUse,iR),yData(validSamples&entriesToUse));    
    
    % Correlation is not valid. Return the data values that were removed
    % and move along
    if isnan(corrVal(iR)) 
        yData = origYData;
        residuals = yData;
        continue
    end
        
    beta(iR,:) = polyfit(xData(validSamples&entriesToUse,iR),yData(validSamples&entriesToUse),modelOrder);
    
    % The polynomial fit is ill formed so ignore and move along
    if any(beta) > maxPolynomialBeta
        beta(iR,:) = nan;
        continue
    end
    
    predict = polyval(beta(iR,:),xData(validSamples,iR));    
    
    yData(validSamples) = yData(validSamples) - predict;
    residuals(:,iR) = yData;

end



