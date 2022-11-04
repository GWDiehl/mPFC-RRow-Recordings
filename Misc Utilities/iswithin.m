function result = iswithin(x,low,hi,inclusive)

% for a matrix of values in x check to see if each is within the low and hi
% boundaries. Low and Hi can be input as vectors in which case it will test
% if each element in x is between any pairing of low-hi entries. 
%
% Variable inclusive is a logical of should the boundary variables be
% included (>=) or not (>).

if ~exist('inclusive','var') || isempty(inclusive)
    inclusive = 1;
end

result = false(size(x));

for w = 1:length(low)
    if inclusive
        r = x>=low(w) & x<=hi(w);
    else
        r = x>low(w) && x<hi(w);
    end
    result = r | result;
end