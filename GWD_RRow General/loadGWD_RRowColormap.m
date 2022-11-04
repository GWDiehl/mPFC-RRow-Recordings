function GWDColorMap = loadGWD_RRowColormap

% Color selections for various mPFC subregions for the GWD mPFC subregions
% manuscript

% GWDiehl 2020


GWDColorMap = struct;


GWDColorMap.All = [0 0 0];

% Original A/B/C grouping
GWDColorMap.GroupA = [249 35 35]/255;
GWDColorMap.GroupB = [165 20 144]/255;
GWDColorMap.GroupC = [32 190 239]/255;

% Various sets of Acc/PL/IL groupings
GWDColorMap.Acc = [232 12 91]/255;
GWDColorMap.PL = [0 137 68]/255;
GWDColorMap.IL = [33 100 219]/255;

GWDColorMap.Acc_TE = GWDColorMap.Acc;
GWDColorMap.IL_TE = GWDColorMap.IL;

% dPL/vPL divisions
GWDColorMap.dPL = [249 145 40]/255;
GWDColorMap.dPL_V2 = GWDColorMap.dPL;
GWDColorMap.dPL_TE = GWDColorMap.dPL;

GWDColorMap.vPL = [47 234 211]/255;
GWDColorMap.vPL_V2 = GWDColorMap.vPL;
GWDColorMap.vPL_TE = GWDColorMap.vPL;

GWDColorMap.No_VO = [.5 .5 .5];
GWDColorMap.GroupVO = [239 143 6]/255;

GWDColorMap.dmPFC = [1 0 0];
GWDColorMap.vmPFC = [0 0 1];

GWDColorMap.Principal = [.8 0 0];
GWDColorMap.Interneuron = [0 0 .8];

GWDColorMap.R506 = 'o';
GWDColorMap.R530 = '+';
GWDColorMap.R535 = '*';
GWDColorMap.R542 = 'x';
GWDColorMap.R536 = '<';
GWDColorMap.R543 = '>';
GWDColorMap.R537 = 'v';
GWDColorMap.R544 = '^';