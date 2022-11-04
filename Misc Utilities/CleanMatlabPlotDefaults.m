function CleanMatlabPlotDefaults

fh1 = figure;
% Sets the axes as black
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{[0,0,0],[0,0,0],[0,0,0]})
set(groot,{'DefaultPolaraxesRColor','DefaultPolaraxesThetaColor'},{[0,0,0],[0,0,0]})
set(groot,'DefaultLegendEdgeColor',[0,0,0])

% This does not seem to work

% % % 
% % % % Sets text color as black
% % % set(groot,'DefaultTextColor',[0,0,0])


close(fh1);