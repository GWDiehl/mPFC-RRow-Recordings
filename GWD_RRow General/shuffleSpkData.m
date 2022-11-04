function shuffledS = shuffleSpkData(S,method,varargin)


% Shuffle the information within spiking data in one of four different
% manners. Timeshift or shuffle ISIs, or if using a Q matrix randomize or
% permute the values. Note, randomize,permute, and timeshift can be used on
% any TSD/cTSD object.
%
% INPUTS
% S - The input data, either a TSD object or a cell array of individual
% cell spike times (TS obj).
% method - Which shuffling method you would like to use.
% 
% Optional
% tStart/tEnd - start/end time of the spiking restriciton for ShuffleISI
% timeBuffer - Buffer on the front/back end (in sec) for circ time shifting
% 
%
% OUTPUT
% shuffledS - Shuffled S data according to the desired method.


% GWD Oct 2020

timeBuffer = 60; % In sec
if isa(S,'tsd') || isa(S,'ctsd')
    tStart = S.starttime;
    tEnd = S.endtime;
else
    tStart = min(cellfun(@starttime, S));
    tEnd = max(cellfun(@endtime, S));
end

commonShiftAcrossCells = 0;

process_varargin(varargin);

if ismember(method,{'Randomize','Permute'}) && ~(isa(S,'tsd') || isa(S,'ctsd'))
    error('Randomize and Permute methods require a tsd/ctsd obj')
elseif strcmp(method,'ShuffISI') && ~strcmp(class(S),'ts')
    error('ShuffleISI method requires raw spike data')
end


switch method
    % Take all of the values and redistribute them. Keeps overall
    % probabilities intact but retains no correspondence across entries.
    case 'Randomize'
        shuffledS = randomizeTSD(S);
        % Jumbles across all dimensions
    case {'Randomize_Full','Randomize_CrossDim'}
        shuffledS = randomizeTSD_CrossDim(S);
        
        % Permute the identity of various entries. Correspondence across
        % identical entries is maintained but overall probabilites are not.
    case 'Permute'
        shuffledS = shuffleTSD(S); 
        % Jumbles across all dimensions
    case {'Permute_Full','Permute_CrossDim'}
        shuffledS = shuffleTSD_CrossDim(S);
        
        % Shuffle the ISI times of spike trains. Only viable for Spike TS
    case 'ShuffISI'
        shuffledS = ShuffleISIs(S, tStart, tEnd);
        
        % Circular timeshift of the data, either raw TS or circshift of a
        % TSD obj. **1 sec resolution on timeshift, uniform distribution**
    case 'TimeShift'
        shuffledS = S;
        shiftAmount = randi(round([timeBuffer,(tEnd - tStart)-timeBuffer]),1);
        
        if isa(S,'tsd') || isa(S,'ctsd')
            shiftAmount = round(shiftAmount/S.dt);
            shuffledS = tsd(shuffledS.range,circshift(shuffledS.data,shiftAmount));
            
        else
            nCells = length(S);               
            if commonShiftAcrossCells
                shiftAmount = repmat(shiftAmount,nCells,1);
            else
                shiftAmount = randi(round([timeBuffer,(tEnd - tStart)-timeBuffer]),nCells,1);
            end
            
            for iC = 1:nCells
                tempTimes = shuffledS{iC}.range + shiftAmount(iC);
                overTime = tempTimes>tEnd;
                tempTimes(overTime) = tempTimes(overTime) - tEnd + tStart;
                shuffledS{iC} = ts(sort(tempTimes));
            end
        end
end

