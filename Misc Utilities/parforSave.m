function parforSave(fileName,variable,modifier)

if nargin<3
    save(fileName,'variable')
else
    save(fileName,modifier,'variable')
end