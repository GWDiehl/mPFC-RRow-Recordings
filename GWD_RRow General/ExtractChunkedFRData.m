function cellFRData = ExtractChunkedFRData(spikingData,triggerTimes,timeEdges,varargin)

% Extract out the binned firing rate data of a set of cell spikings around
% a set of trigger times. This is largly equivalent to making a PETH, but
% allows for uneven timeEdges, and retains the avg FR for each cell in each
% bin on each event as opposed to just the avg across all events. Does this
% all by building a Q Mtx around the timing bins of interest.
%
% One useful mechanism is that this code allows for building around unequal
% time edges either unequal width bins for a given event or uneven across
% events

% GWD May 2022


% How do we want to chunk the event times for grouping? Every X events.
% Make sure there is no overlap between time windows so use larger
% chunkSteps if there is likelyhood for windows to overlap between events.
chunkStep = 5;

itterCount = 1; % Check to make sure that the fuction cannt recursive itself to inf

process_varargin(varargin);

if chunkStep == 1
    chunkStep = false;
end

nCells = length(spikingData);
nEvents = length(triggerTimes);
nTimes = size(timeEdges,2)-1;

% Grand output FR (Cells x Time x Events)
cellFRData = nan(nCells,nTimes,nEvents);

if nEvents == 0 || nTimes == 0
    return
end

if size(timeEdges,1)==1
    timeEdges = repmat(timeEdges,nEvents,1);
end
timeCenters = cell2mat(arrayfun(@(x) timeEdges(x,2:end) - diff(timeEdges(x,:))/2,[1:nEvents]','UniformOutput',0));


tStart = triggerTimes(1)+timeEdges(1)-diff(timeEdges(1,1:2));
tEnd = triggerTimes(end)+timeEdges(end)+diff(timeEdges(end,end-1:end));



if chunkStep
    % One option whereby event times are grouped together evey
    % several events to require fewer builds of a Q mtx. Care needs
    % to be taken however to ensure that the chunk step size is
    % large enough that event windows do not overlap. Dramatically
    % reduces processing time.
    
    chunkedEvents = nan(chunkStep,ceil(nEvents/chunkStep));
    chunkedEvents(1:nEvents) = 1:nEvents;
    
    for iE = 1:size(chunkedEvents,1)
        tempIdx = chunkedEvents(iE,:);
        tempIdx(isnan(tempIdx)) = [];
        if isempty(tempIdx)
            continue
        end
        % Collect up all of the requisit QTimes (buffer one extra
        % bin on either side for safty)
        QTimes = [triggerTimes(tempIdx)+...
            [timeEdges(tempIdx,1)-diff(timeEdges(tempIdx,1:2),[],2),...
            timeEdges(tempIdx,:),...
            timeEdges(tempIdx,end)+diff(timeEdges(tempIdx,end-1:end),[],2)]]';
        
        if ~issorted(QTimes(:))
            warning('You have overlap in your timing, retrying with a larger chunk size'); 
            if itterCount > 10
                error('You have gone through 5 itterations and still having problems, something is wrong, fix it')
            end
            cellFRData = ExtractChunkedFRData(spikingData,triggerTimes,timeEdges,'chunkStep',chunkStep*2,'itterCount',itterCount+1);
            return
        end
        QTemp = MakeQfromS_NonUniformTime(spikingData, ts(QTimes(:)),'tStart',tStart,'tEnd',tEnd);
        QTemp = recenterQMtx(QTemp);
        for iF = 1:length(tempIdx)
            triggIdx = tempIdx(iF);
            currEventTimes = triggerTimes(triggIdx) + timeCenters(triggIdx,:);
            cellFRData(:,:,triggIdx) = full(QTemp.data(currEventTimes))'./diff(timeEdges(triggIdx,:));
        end
    end
else
    % Deal with each event entriley independently
    for iE = 1:nEvents
        currEventTimes = triggerTimes(iE) + timeCenters(iE,:);
        
        QTimes = triggerTimes(iE)+[timeEdges(iE,1)-diff(timeEdges(iE,1:2)),timeEdges(iE,:),timeEdges(iE,end)+diff(timeEdges(iE,end-1:end))];
        QTemp = MakeQfromS_NonUniformTime(spikingData, ts(QTimes),'tStart',sessStart,'tEnd',sessEnd);
        QTemp = recenterQMtx(QTemp);
        
        % Reshape the Shuff back to D2
        cellFRData(:,:,iE) = full(QTemp.data(currEventTimes))'./diff(timeEdges(iE,:));
    end
end

