function [fh, thresholds] = plotRRowChoices(delay,accept,quit,varargin)

% Plot a figure of the RRow choices split by restraurnat and identifying
% accept/skip/quit as a function of delay. Optional plot a threshold value
% (default: Heaviside; Options: None, Sigmoid)
%
% Inputs
% delay - A Laps x Sites matrix of the delay at each site (Nans OK)
% accept - A logical array of all laps that were accepted (NO Nans)
% quit - A logical array of all laps that were quit
%
% Output
% fh - Figure handle of the resultant plot
% thresholds - A 1 x Sites vector of the computed threhsolds

% GWD Jan 2021!!!

%%

thresholdType = 'Heaviside';
process_varargin(varargin);

%%

earned = accept & ~quit;

accept = double(accept);
accept(quit) = nan;

minDelay = min(delay(:));
maxDelay = max(delay(:));


nSites = size(delay,2);
if strcmp(thresholdType,'Both')
    thresholds = nan(2,nSites);
else
    thresholds = nan(1,nSites);
end

fh = figure;
[sub1, sub2] = computeSubplotSize(nSites);

for iZ = 1:nSites
    validLaps = ~isnan(delay(:,iZ));
    
    curDelay = delay(validLaps,iZ);
    curQuit = quit(validLaps,iZ);
    curAccept = accept(validLaps,iZ);
    curEarn = earned(validLaps,iZ);
    
    subplot(sub1,sub2,iZ)
    hold on
    plot(curDelay+rand(size(curDelay))/4-.125,curAccept+rand(size(curAccept))/4-.125,'o')    
    plot(curDelay(curQuit)+rand(sum(curQuit))/4-.125,.5+rand(sum(curQuit))/4-.125,'.k','MarkerSize',20)
    
    axis([minDelay-3 maxDelay+3 -.25 1.25])
    
    switch thresholdType
        case 'Heaviside'
            thresholds(iZ) = RRheaviside(curDelay,curEarn);
%             plotVertLine(thresholds(iZ),{'k--'})
            plot([minDelay,thresholds(iZ),thresholds(iZ),maxDelay],[1,1,0,0],'k')
            
        case 'Sigmoid'
            [param,stat]=sigm_fit(curDelay,curEarn,[0 1 nan nan],[nan nan nanmedian(curDelay) -1],0);
            dataIn = linspace(minDelay,maxDelay,100);
            dataOut = sigmaFunction(dataIn,param);
            plot(dataIn,dataOut,'r--')
            plot(param(3),.5,'m*')
            thresholds(iZ) = param(3);
        case 'Both'
            thresholds(1,iZ) = RRheaviside(curDelay,curEarn);
            %             plotVertLine(thresholds(1,iZ),{'k--'})
            plot([minDelay,thresholds(1,iZ),thresholds(1,iZ),maxDelay],[1,1,0,0],'k')
            
            [param,stat]=sigm_fit(curDelay,curEarn,[0 1 nan nan],[nan nan nanmedian(curDelay) -1],0);
            dataIn = linspace(minDelay,maxDelay,100);
            dataOut = sigmaFunction(dataIn,param);
            plot(dataIn,dataOut,'r--')
            plot(param(3),.5,'m*')
            thresholds(2,iZ) = param(3);
        case 'None'
            % All good
        otherwise
            error('Unknown threshold type')
    end   
    
    title(['Restaurant ',num2str(iZ)])
    xlabel('Delay')
    ylabel('Earn')
end
    
    