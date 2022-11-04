function RecenteredQ = recenterQMtx(Q)

% Takes an original Q matrix where time values represent the beginning of
% the time bin and recenters the data such that times represent the mid
% point of the time bin. Note this will reduce the time dimention by 1 bin.

% GWD 2020/9/15


dataType = class(Q);

oldData = Q.data;
oldTime = Q.range;

newData = oldData(1:end-1,:);
newTime = oldTime(1:end-1)+diff(oldTime)/2;

switch dataType
    case 'tsd'
        RecenteredQ = tsd(newTime,newData);
    case 'ctsd'
        RecenteredQ = ctsd(newTime(1),Q.dt,newData);
end