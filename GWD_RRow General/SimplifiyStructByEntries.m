function simplifiedStruct = SimplifiyStructByEntries(struct,selectedField,varargin)

% For a multilevel structure, select out a single data field and make that
% the only entry of the immediate parent strucutre

process_varargin(varargin);

values = fieldnames(struct);
for v = 1:length(values)
    if ~ismember(selectedField,fieldnames(struct.(values{v})))
        simplifiedStruct.(values{v}) = SimplifiyStructByEntries(struct.(values{v}),selectedField);
    else
        simplifiedStruct.(values{v}) = struct.(values{v}).(selectedField);
    end
end