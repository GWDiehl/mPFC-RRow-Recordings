function [r, p] = nancorr(A,B)

vectorA = reshape(A,[],1);
vectorB = reshape(B,[],1);

usable = ~isnan(vectorA) & ~isnan(vectorB);

if ~any(usable)
    r = nan;
    p = nan;
    return
end

[r, p] = corr(vectorA(usable),vectorB(usable));


