function groupedData = groupRRowBySite(R,varargin)

sitesOfInterest = unique(R.ZoneIn);
fieldsOfInterest = {'EnteringZoneTime' 'ExitZoneTime' 'ZoneDelay' 'SkipOffer' 'QuitOffer' 'CurrentCycle' 'CurrentLap' 'FeederTimes'};

maxLaps = 125;

process_varargin(varargin);

nSites = length(sitesOfInterest);

groupedData = struct;

for f = 1:length(fieldsOfInterest)
    groupedData.(fieldsOfInterest{f}) = nan(maxLaps,nSites);
    
    switch fieldsOfInterest{f}
        case {'EnteringZoneTime' 'ExitZoneTime' 'ZoneDelay' 'SkipOffer' 'QuitOffer' 'CurrentCycle' 'CurrentLap'}
            for z = 1:length(sitesOfInterest)
                toUse = R.ZoneIn == sitesOfInterest(z);
                groupedData.(fieldsOfInterest{f})(R.CurrentCycle(toUse),z) = R.(fieldsOfInterest{f})(toUse);
            end
        case {'FeederTimes'}
            for z = 1:length(sitesOfInterest)                
                % Populate Feeder times into to full lap matrix
                FeederTimes = R.FeederTimes(R.FeedersFired == sitesOfInterest(z));
                toUse = R.ZoneIn == sitesOfInterest(z);                
                earnLap = R.QuitOffer(toUse) == 0;
                lapNumber = R.CurrentCycle(toUse);
                groupedData.FeederTimes(lapNumber(earnLap),z) = FeederTimes;
            end
    end
end
    
    
    
