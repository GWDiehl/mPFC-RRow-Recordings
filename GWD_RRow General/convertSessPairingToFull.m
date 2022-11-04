function pairingFull = convertSessPairingToFull(pairingSess)

% Take a series of by session pairing of entries and group them all
% together into a single comprehensive set. Requires adjusting the
% respective numbers accordingly.
%
% Input
% pairingSess - Cell array of a N x D (commonly nCells x 2) mtx of pairwise
%       relationships
%
% Output
% pairingFull - N x D matrix where N is the total number of pairs across
%       all sessions adjusted for each per session count

% GWD Oct 2020


% Adapt the pairing numbers so they now match a full cat of the data
counts = cellfun(@(x) length(unique(x)),pairingSess);
runningCounts = [0; cumsum(counts(1:end-1))];

% Pairings now reflect adjusted position in the full cat list
pairingFull = cellfun(@(x,y) x+y, pairingSess,num2cell(runningCounts),'UniformOutput',0);

pairingFull = cell2mat(pairingFull);
