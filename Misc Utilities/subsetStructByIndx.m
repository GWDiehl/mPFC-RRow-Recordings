function subsetStruct = subsetStructByIndx(struct,indx,varargin)

% For a structure composed of a series of matrices, take a subset of each
% field matrix based on an index selected along the Mtx first dimension. In
% the case of only a sinlge indx selected you can squeeze out the first
% dimension based on a varargin input.

% GWD 2019

squeezeD1 = 0;

process_varargin(varargin);

if islogical(indx)
    indx = FastFind(indx);
end

values = fields(struct);
for v = 1:length(values)
    if isstruct(struct.(values{v}))
        subsetStruct.(values{v}) = subsetStructByIndx(struct.(values{v}),indx);
    else
        subsetStruct.(values{v}) = selectalongfirstdimension(struct.(values{v}),indx);
        
        % If we dealing down to a single indx we can compress out the first
        % dimention
        if length(indx) == 1 && squeezeD1
            subsetStruct.(values{v}) = squeeze(subsetStruct.(values{v}));
        end
    end
end