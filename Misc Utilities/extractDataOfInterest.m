function [data] = extractDataOfInterest(data,indx)


dataFields = fields(data);
for f = 1:length(dataFields)
    if istsd(data.(dataFields{f}))
        data.(dataFields{f}).D = data.(dataFields{f}).D(:,indx);
    else
        data.(dataFields{f}) = data.(dataFields{f})(indx);
    end
end