function [Threshold,Slope] = RRowSigmoidFit(Delay,StayGo,varargin)


boundMax = 1;
boundMin = 0;
boundMedian = nan;
boundSlope = nan;

initialMax = 1;
initialMin = 0;
initialMedian = nanmedian(Delay);
initialSlope = -1;

process_varargin(varargin);

usableData = ~isnan(Delay) & ~isnan(StayGo); % Dont include nans

% There must be more than 1 delay to calculate a threshold
if length(unique(Delay(usableData)))<=1 
    Threshold = nan;
    Slope = nan;
    return
end


boundConditions = [boundMin,boundMax,boundMedian,boundSlope];
initialConditions = [initialMin,initialMax,initialMedian,initialSlope];
[fitParams,stat]=sigm_fit(Delay(usableData),StayGo(usableData),boundConditions,initialConditions,0);

fitMin = fitParams(1);
fitMax = fitParams(2);
Threshold = fitParams(3);
Slope = fitParams(4);

if Threshold > max(Delay) || Threshold < min(Delay) % Outside of range so bad fit.
    Threshold = nan;
    Slope = nan;
end

