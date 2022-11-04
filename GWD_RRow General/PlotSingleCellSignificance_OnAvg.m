function [inputFig] = PlotSingleCellSignificance_OnAvg(inputFig,sigTime,dataCenters,anatGroup,sigYPos,varargin)

GWDColor = loadGWD_RRowColormap;

process_varargin(varargin);

nAnat = length(anatGroup);

nPlots = size(sigTime,3);
[subA,subB] = computeSubplotSize(nPlots);

figure(inputFig)
for iP = 1:nPlots
    subplot(subA,subB,iP)
    hold on
    for iA = 1:nAnat
        plot(dataCenters,sigTime(iA,:,iP)*sigYPos(iA,iP),'o','color',GWDColor.(anatGroup{iA}))
    end
end

