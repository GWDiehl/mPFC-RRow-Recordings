function SelectDefaultColorMap(colormap)

tempFig = figure;
set(groot,'DefaultFigureColormap',colormap);
close(tempFig);