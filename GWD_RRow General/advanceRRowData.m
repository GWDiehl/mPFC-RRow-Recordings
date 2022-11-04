function advancedData = advanceRRowData(data,movementSpots,varargin)

% Advances RRow mtx data through to the next/previous site. Good for
% finding a particular value at the subsequent(-)/preceeding(+) site
% reletive to current.
%
% INPUTS 
% data - A 2D/3D mtx of RRow data organized as Sess x Laps x Sites. If
%   using 2D it is expected that the Sess dim has been removed (only one
%   sess)
% movementSpots - How many places forwards(+) or backwards(-) you want to
%   advance the input data. Use Postive values to look at the previous
%   visit's data and negative values to look at the subsequent visit's data
%
% Optional:
% retainSite - When advancing do you want to maintain data within a site.
%   I.E. data will advance to the next entries at the same RR site. (default
%   false, advance to next sequental site regardless of RR)
%
% Outputs 
% advancedData - Mtx of the same size as input data but now with entries
%   moved accroding to movementSpots

% GWD 2019

retainSite = false;

process_varargin(varargin);

% If the data are only 2D there is no Sess Dim, only Laps x Site
nDim =  ndims(data);
if nDim == 2
    data = reshape(data,1,size(data,1),size(data,2));
end

switch retainSite
    case true
        advancedData = circshift(data,movementSpots,2);
    case false
        [nSess, nLaps, nSites] = size(data);
        
        advancedData = data;
        
        for s = 1:nSess
            % Extract out session data and transpose so it is now Sites x Laps
            sessionData = squeeze(data(s,:,:))';
            % Vecotrize it
            vectorData = sessionData(:);
            % Advance it
            shiftedData = circshift(vectorData,movementSpots);
            % Reformat back into Sites x Laps
            reformatedData = reshape(shiftedData,nSites,nLaps);
            % Transpose back to the original and populate the output matrix
            advancedData(s,:,:) = reformatedData';
        end
end

% Return to 2D if no Sess field
if nDim == 2
    advancedData = squeeze(advancedData);
end
