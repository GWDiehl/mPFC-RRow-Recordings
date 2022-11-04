function [basicQuant,defaultRange] = QuantifyBasicCellSpiking(spikes,selectedQuant,cellWF,varargin)


% GWD Nov 2021

sessDuration = 60*60;
Q_dt = .05;

shortISI = .01;
longISI = .125;

WF_SampleRate = 30000;
WFPeak = [4 12];

process_varargin(varargin);

%%

nSpks = cellfun(@(x) length(x.range),spikes);
nCells = length(spikes);

if ismember(selectedQuant,{'SpikeWidth' 'PeakValleyRatio'})
    assert(~isempty(cellWF),'You have not provided the WF data')
end

basicQuant = nan(length(spikes),1);    

switch selectedQuant
    case 'AvgRate'        
        basicQuant = nSpks/sessDuration;
        defaultRange = [-.5 40];
        
    case 'MedianISI'
        basicQuant = cellfun(@(x) x.dt,spikes);
        defaultRange = [-.05 2.5];
        
    case 'ISIStd'
        basicQuant = cellfun(@(x) nanstd(diff(x.range)),spikes);
        defaultRange = [-.5 20];
        
    case {'ISIRatio' 'ISIRatio_Log'}
        shortISIPort = cellfun(@(x) sum(diff(x.range) < shortISI),spikes)./nSpks;
        longISIPort = cellfun(@(x) sum(diff(x.range) > longISI),spikes)./nSpks;
        
        basicQuant = shortISIPort./longISIPort;
        defaultRange = [-.05 5];
        
        if strcmp(selectedQuant,'ISIRatio_Log')
            basicQuant = log(basicQuant);
            defaultRange = [-8 5];
        end
        
    case 'COV'
        QMtx = MakeQfromS(spikes,Q_dt);
        spikingData = QMtx.data/Q_dt;
        for iC = 1:nCells
            basicQuant(iC) = nanstd(spikingData(:,iC))/nanmean(spikingData(:,iC));
        end
        defaultRange = [0 50];
        
    case 'SpikeWidth'
        for iC = 1:nCells
            tempWF = cellWF{iC};
            if ~any(~isnan(tempWF))
                continue
            end
            tempWF(1:WFPeak(1)) = nan;
            [~, peakIdx] = min(tempWF);
            tempWF(1:peakIdx) = nan;
            [~,vallyIdx] = max(tempWF);
            
            % Convert from samples to ms
            basicQuant(iC) = (vallyIdx-peakIdx)/WF_SampleRate*1000;
        end        
        defaultRange = [0 1];
        
    case 'PeakValleyRatio'
        for iC = 1:nCells
            tempWF = cellWF{iC};
            if ~any(~isnan(tempWF))
                continue
            end
            tempWF(1:WFPeak(1)) = nan;
            [peak, peakIdx] = min(tempWF);
            tempWF(1:peakIdx) = nan;
            vally = max(tempWF);
            
            basicQuant(iC) = -peak/vally;
        end        
        defaultRange = [1 9];
        
    otherwise
        error('Case not yet implimented')
end

