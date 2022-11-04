function [IdPhi, AvgdPhi] = compute_IdPhi_AvgdPhi(x, y, tstart, tend, varargin)

% IdPhi = compute_IdPhi(x, y, tstart, tend, varargin)
%
% Calculates zIdPhi
% 2012-04-16 AndyP added arbitrary (tstart, tend) times
% 2012-06-21 AndyP fixed x,y assignment error line 17
% 2012-06-21 AndyP added VT1 or VT2 option
%
% 2012-07-25 ADR tend defined as ExitingCPTime
% 2013-01-21 AndyP added checks for sorted timestamps, dxdt returns error if timestamps are out of order  

% 2018-11-21 GWD Extracted the computations from zIdPhi to receive any x,y
% and CP times outside of the standard sd format.

% 2019-01-18 GWD Added an AvgdPhi field which will calculate the avg change
% in phi metric over the period of interest as opposed to integrating the
% signal; This is important for cases in which the total time/distance over
% which dPhi is calcualted is not consistant across conditions/comparisons
% (For example different distances covered for different choices or different time
% to make them)

dxdtWindow = 0.5;
dxdtSmoothing = 0.1;

process_varargin(varargin);
assert(~isempty(tstart), 'Unknown CP entry times.');
assert(~isempty(tend), 'Unknown CP exit times.');

%%%%%%%%%%%%%%%%%%%%
assert(length(tstart)==length(tend),'tstart must equal tend');
if ~issorted(x.range); time=sort(x.range); x=tsd(time,x.data); end % cheetah (ring buffer error?) causes out of order timestamps
if ~issorted(y.range); time=sort(y.range); y=tsd(time,y.data); end
[ dx ] = dxdt(x, 'window', dxdtWindow, 'postSmoothing',dxdtSmoothing);
[ dy ] = dxdt(y, 'window', dxdtWindow, 'postSmoothing',dxdtSmoothing);

phi = tsd(dx.range(), atan2(dy.data(), dx.data()));
uphi = tsd(phi.range(), unwrap(phi.data()));
dphi = dxdt(uphi, 'window', dxdtWindow, 'postSmoothing',dxdtSmoothing);
%%%%%%%%%%%%%%%%%%%%
nPasses = length(tstart);
IdPhi = nan(nPasses,1);
AvgdPhi = nan(nPasses,1);

for iL = 1:nPasses
    if isnan(tstart(iL))
        continue
    end
	dphi0 = dphi.restrict(tstart(iL), tend(iL));
	IdPhi(iL) = sum(abs(dphi0.data()));
    AvgdPhi(iL) = nanmean(abs(dphi0.data()));
end




