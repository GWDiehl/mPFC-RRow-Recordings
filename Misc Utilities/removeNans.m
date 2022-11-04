function newdata = removeNans(data,dim)

if ~exist('dim','var') || isempty(dim)
    dim = 1;
end

switch dim
    case 1
        newdata = data(:,~isnan(nanmean(data)));
    case 2
        newdata = data(~isnan(nanmean(data,2)),:);
end
