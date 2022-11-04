function [tsd_out, mask_out] = restrictdata_CrossTSD(tsa_in, tsd_restrictor, dmin, dmax, varargin)

% Takes an input TSD and restricts its T/D pairings to epochs in which a
% second TSD falls within a desired data range. Allows for restricting
% TSD_A based on the data in TSD_B.
%
% Inputs
% tsa_in - The primary TSD/TS object that will have its data
%       restricted/masked
% tsd_restrictor - The TSD whose data will be used for the purpose of
%       restriting
% dmin/dmax - The data min and max for which we are considering the entries
%       in tsd_restrictor. Of interest are the epochs during with tsd-restrictor
%       falls within this range.
%
% Optional Inputs
% maskData - Do we want to simply mask out the unwanted data. Default -
%       True. Alternativly an entry of false will cause the wanted data to
%       be completely disregarded from the tsd.
% excludeWindow - Instead of keeping only those epochs within the desired
%       data range we will disregard these epochs. Default false (keep the
%       data within the window)
%
% Outputs
% tsd_out - The desired TSD/TS in which the tsa_in has been
%       restricted/masked
% mask_out - A TSD object of the mask applied to the original tsa_in

% GWD April 2020

maskData = true;
excludeWindow = false;

process_varargin(varargin);

assert(isa(tsa_in, 'ts')|| isa(tsa_in, 'ctsd'),'Data must be a TS/TSD object')

% Identify the valid time entries in the restrictOR TSD
[~, validTimes] = restrictdata(tsd_restrictor,dmin,dmax);

% For each time point in tsa_in see if it was within the valid window
% TSD.data does not handle logicals well so requires extrapolating to
% maintain sizing.
mask = data(validTimes,tsa_in.range,'extrapolate',1);

% If we want to exclude pts in the window (keep those out of window) flip
% the mask
if excludeWindow
    mask = ~mask;
end

% Any mask entries derived from extrapolation should be excluded from the
% mask as these were outside of the relevant data to begin with.
% % % (restrictEE data points fall outside of the times covered by the
% % % restritOR)
windowStart = tsd_restrictor.starttime;
windowEnd = tsd_restrictor.endtime;
inWindow = tsa_in.range >= windowStart & tsa_in.range <=windowEnd;
mask(~inWindow) = 0;


mask_out = tsd(tsa_in.range,mask);

% For TSD, nan out those pts outside of the mask
if isa(tsa_in,'tsd') || isa(tsa_in,'ctsd')
    if maskData
        D = tsa_in.data;
        D(~mask,:) = nan;
        
        tsd_out = tsd(tsa_in.range,D);
        
        % Or disregard them from the TSD entriely
    else
        D = selectalongfirstdimension(tsa_in.data,mask);
        T = tsa_in.range;
        tsd_out = tsd(T(mask),D);
    end
    
    % For TS objects only keep those entries within the mask
else
    T = tsa_in.range;
    tsd_out = ts(T(mask));
end

