function dataAreUnique = isunique(testData);

% Check to see if the data are unique. Get the unique elements and compare
% the length to the original. All Nan entries will be treated independnetly
% as default behavior for 'unique'.

% GWD June 2022

uniqueData = unique(testData);
dataAreUnique = isequal(length(uniqueData),length(testData));

