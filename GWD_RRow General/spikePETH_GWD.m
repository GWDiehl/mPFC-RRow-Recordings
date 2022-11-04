function [PETH, PETH_Var] = spikePETH_GWD(S,t,varargin)

% spikePETH(S, t, varargin)
%
% input:
%        S - cell array of spike TS objects
%        t - Vector of event times
% Optional:
%       window - Window size in sec
%       dt - Bin width in sec
%           bounds - Additional restriction bounds around each t for limiting
%       PETH window. High utility if restricting data but where overlapping
%           PETH windows limit restriction efficecy.
%
% output:
%        TSD PETH - tsd with time as window bin centers, data as sum in
%           each bin
%        TSD PETH_Var - tsd with time as window bin centers, std of the sum
%           across events
%
% parameters:
%   window = [-2 5]
%   dt = 0.01; % seconds
%

% GWD April 2020 Adapted from spikePETH to simplify inerworkings/outputs,
% handle a cell array of spiking data, and compute the var(Std
% specifically) of counts across events

% NOTE: This version outputs times as bin centers with nCenters = nEdges-1
% whereas the original spikePETH code output time as left bin edge with the
% very last data matching only the edge time not a full bin width.

%--------------------------
% parameters
%--------------------------
window = [-2 5];
dt = 0.01;
bounds = [-inf inf];

process_varargin(varargin);

%-------------------------
% prep
%--------------------------
nT = length(t);
nC = length(S);

if size(bounds,1) == 1
    bounds = repmat(bounds,nT,1);
end

x = window(1):dt:window(2);
binCenters = x(1:end-1) + diff(x)/2;

PETH = cell(nC,1);
PETH_Var = cell(nC,1);


%---------------------------
% go
%---------------------------


for iC = 1:nC
    outputS = [];
    
    V0 = nan(length(x),nT);
    
    % Calculate the SpikeTimes as deviation from EventTime
    for iT = 1:nT
        startBound = max(window(1),bounds(iT,1));
        endBound = min(window(2),bounds(iT,2));
        
        S0 = data(S{iC}.restrict(t(iT)+startBound, t(iT)+endBound));
        S0 =  S0-t(iT);
        outputS = cat(1, outputS, S0);        
        
        % Empty input to histc returns empty, not zeros
        V0(:,iT) = histc([nan;S0], x); % Count each 
    end
    
    % Compute Full Count and Std of per event count
    % Empty input to histc returns empty, not zeros
    H = histc([nan;outputS], x);    
    V = nanstd(V0,[],2);
    
    % Make sure it is a column vector
    if size(H,2)>1
        H = reshape(H,[],1);
    end
    
    % Very last bin of H is count == final bin edge (quirk in histc)    
    PETH{iC} = tsd(binCenters', H(1:end-1));
    PETH_Var{iC} = tsd(binCenters', V(1:end-1));    
end







