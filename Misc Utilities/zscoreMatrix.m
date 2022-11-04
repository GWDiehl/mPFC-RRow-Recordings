function newMatrix = zscoreMatrix(matrix,dim)

if ~exist('dim','var') || isempty(dim)
    dim = 1;
end

newMatrix = matrix;

elements = size(newMatrix,dim);
switch dim
    case 1        
        for e = 1:elements
            newMatrix(e,:) = (newMatrix(e,:)-nanmean(newMatrix(e,:)))/nanstd(newMatrix(e,:));
        end        
    case 2
        for e = 1:elements
            newMatrix(:,e) = (newMatrix(:,e)-nanmean(newMatrix(:,e)))/nanstd(newMatrix(:,e));
        end
end