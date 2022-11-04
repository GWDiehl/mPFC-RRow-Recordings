function subDataDef = subsetFromDataDef(dataDef,Sess)

% Extracts a subset of sessions from the dataDef strucure

% GWD 2019

dataDefFields = fields(dataDef);

switch length(Sess)
    
    % If you only want one session then pull out things that are in cells
    case 1
        for f = 1:length(dataDefFields)
            if iscell(dataDef.(dataDefFields{f}))
                subDataDef.(dataDefFields{f}) = dataDef.(dataDefFields{f}){Sess};
            else
                subDataDef.(dataDefFields{f}) = dataDef.(dataDefFields{f})(Sess);
            end            
        end
    % Otherwise keep the dataDef as is but only hold onto a subset of it  
    otherwise
        for f = 1:length(dataDefFields)
            subDataDef.(dataDefFields{f}) = dataDef.(dataDefFields{f})(Sess);
        end
end