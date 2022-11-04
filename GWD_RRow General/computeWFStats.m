function [WFStats] = computeWFStats(WF,varargin)

[~, channel] = max(range(WF));
[peak, riseTime] = nanmax(WF(:,channel));
[valley, minTime] = nanmin(WF(:,channel));
if peak < valley
    WF = -WF;
    [peak, riseTime] = nanmax(WF(:,channel));
    [valley, minTime] = nanmin(WF(:,channel));
end
WF = WF(:,channel);


fallTime = minTime - riseTime;
riseTime = riseTime*(1000/32); %convert from bins to microseconds based on our 32k sampling rate
fallTime = fallTime*(1000/32);

ratios = abs(peak/valley);


% Code written my E Mankin (exercise) that I cant understand to be able to
% rewrite and simplify
iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
paren = @(x, varargin) x(varargin{:});
posSignchange = @(x)find(diff(sign(x))>0);
negSignchange = @(x)find(diff(sign(x))<0);
minima = @(x)posSignchange(diff([max(x)+1,x,max(x)+1]));
maxima = @(x)negSignchange(diff([min(x)-1,x,min(x)-1]));
absMinHelper = @(x)find(x(minima(x))==min(x));
absMaxHelper = @(x)find(x(maxima(x))==max(x));
wavewidth3 = @(x)iif(isempty(minima(x))||isempty(maxima(x)),@()NaN,...
    true,@()(paren(minima(x),absMinHelper(x))-paren(maxima(x),absMaxHelper(x)))/32);


width3 = cellfun(@(x)wavewidth3(x),{WF'});


WFStats.WF = WF;
WFStats.width3 = width3;
WFStats.ratios = ratios;
WFStats.peak = peak;
WFStats.valley = valley;
WFStats.riseTime = riseTime;
WFStats.fallTime = fallTime;
    
    