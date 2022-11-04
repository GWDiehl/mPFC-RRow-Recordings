function mtx = NanMatrixDiag(mtx)

% Replaced the diagonals of a Matrix with NaNs

mtx(1:size(mtx,1)+1:end) = nan;