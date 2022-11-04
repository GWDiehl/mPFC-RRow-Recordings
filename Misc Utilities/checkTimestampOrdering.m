function [NewPathData, NewX, NewY] = checkTimestampOrdering(PathData,x,y)

% Find if timestamps in a set of tracking data are out of order. In the
% case that there are values that preceed the t-1 entry, simply remove
% those data. In the case that these instances represent more than 5% of
% the total, error out because you have something serioulsy wrong.

% GWD 2019

if nargin == 3
    t = PathData;
    
    outOfOrder = [0; diff(t)<0];
    if outOfOrder(2)
        outOfOrder(1) = 1;
        outOfOrder(2) = 0;
    end
    
    if sum(outOfOrder)/length(outOfOrder) > .05
        error('More than 5% of samples are out of order, you have a bigger problem')
    end
    
    NewPathData = t(~outOfOrder);
    NewX = x(~outOfOrder);
    NewY = y(~outOfOrder);
    
    % This would be the more complicated situation of not just passing in
    % t,x,y values but rather some fancy structured configuration of
    % things.
else
    
    switch class(PathData)
        case 'tsd'
            % Which timepoints are out of order
            outOfOrder = [0; diff(PathData.T)<0];
            
            % If the second timestamp is out of order chances are the first sample is the error
            if outOfOrder(2)
                outOfOrder(1) = 1;
                outOfOrder(2) = 0;
            end
            if sum(outOfOrder)/length(outOfOrder) > .05
                error('More than 5% of samples are out of order, you have a bigger problem')
            end
            
            NewPathData = tsd(PathData.T(~outOfOrder),PathData.D(~outOfOrder));
            
        case 'struct'
            variables = fields(PathData);
            NewPathData = struct;
            
            if isa(PathData.(variables{1}),'tsd') % If it is a collection of TSDs
                for i = 1:length(variables)
                    tempData = PathData.(variables{i});
                    
                    outOfOrder = [0; diff(tempData.T)<0];
                    % If the second timestamp is out of order chances are the first sample is the error
                    if outOfOrder(2)
                        outOfOrder(1) = 1;
                        outOfOrder(2) = 0;
                    end
                    if sum(outOfOrder)/length(outOfOrder) > .05
                        error('More than 5% of samples are out of order, you have a bigger problem')
                    end
                    
                    NewPathData.(variables{i}) = tsd(tempData.T(~outOfOrder),tempData.D(~outOfOrder));
                end
                
            else % Otherwise it is presumibly just vectors of values
                if ismember(variables,'t')
                    outOfOrder = [0; diff(PathData.t)<0];
                elseif ismember(variables,'T')
                    outOfOrder = [0; diff(PathData.T)<0];
                else
                    error('How are timestamps identified??')
                end
                
                if outOfOrder(2)
                    outOfOrder(1) = 1;
                    outOfOrder(2) = 0;
                end
                
                if sum(outOfOrder)/length(outOfOrder) > .05
                    error('More than 5% of samples are out of order, you have a bigger problem')
                end
                
                for i = 1:length(variables)
                    NewPathData.(variables{i}) = PathData.(variables{i})(~outOfOrder);
                end
            end
            
            
    end
end
