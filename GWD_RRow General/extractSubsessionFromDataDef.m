function subDataDef = extractSubsessionFromDataDef(dataDef,Sess)



dataDefFields = fields(dataDef);
subDataDef = struct;

switch length(Sess)
    
    case 1
        for f = 1:length(dataDefFields)
            if iscell(dataDef.(dataDefFields{f}))
                subDataDef.(dataDefFields{f}) = dataDef.(dataDefFields{f}){Sess};
            else
                subDataDef.(dataDefFields{f}) = dataDef.(dataDefFields{f})(Sess);
            end            
        end
        
    otherwise
        for f = 1:length(dataDefFields)
            for s = 1:length(Sess)
                if iscell(dataDef.(dataDefFields{f}))
                    subDataDef.(dataDefFields{f}){s} = dataDef.(dataDefFields{f}){s};
                else
                    subDataDef.(dataDefFields{f})(s) = dataDef.(dataDefFields{f})(s);
                end
            end
        end
end