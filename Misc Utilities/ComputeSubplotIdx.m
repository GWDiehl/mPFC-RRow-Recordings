function subplotIdx = ComputeSubplotIdx(subA,subB,idxA,idxB)

% Computes the counted index value of a plot coordinate in subplots
% GWD April 2021

if idxA > subA || idxB > subB
    error('Your requested position does not exist in this subplot dimensions')
end

subplotIdx = (idxA-1)*subB + idxB;