function [Location, Depth, Chan, Dist, SelectedFNs,fh1] = LocalizeUnits_mPFC(varargin)

% Standalone function that will go through and load all of the mPFC RRow
% (standard Task) data from GWD project. Relies on 'FindUnitLocation_mPFC'
% and basically just serves as a useful means of gathering locations of ALL
% cells independent of anything else.
%
% Optional Input:
% plotData - Plots a 3D scatter plot of all of the recorded cells, colored
% by rat. (Default on).
%
% Outputs:
% Location - An nCell x 3 Matrix of the locaiton of each unit in
%       standard space.
% Depth - The DV distance in real space that each unit is located from
%       the surrace of brain
% Chan - The sorted channel where each unit was localized on its
%       respective probe.
% Dist - The distance in real space of the unit from the bottom of the
%       probe.
% Full_FNs - The file name of each respective unit included

% GWD 2020

% Directories
[rootDir, DataDefName] = ID_GWD_Dirs;
dataDef = load([rootDir,DataDefName]);
DataDefOI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

% Grab the sub directories
subDataDef = subsetFromDataDef(dataDef,DataDefOI);
SSNs = subDataDef.SSN;

plotData = 1;
axisVals = [-5 5 2 4.5 -7 0];

anatGroup = 'All';
anatplotType = 'Sections';

process_varargin(varargin);

%%
% Data Def loading
nSess = length(SSNs);

Full_FNs = cell(nSess,1);

%% For each session loop through and do stuff

for iS = 1:nSess
    [~,~,~, Full_FNs{iS}] = GWD_LoadRRowSession(SSNs{iS},{'SpikeNames'});
end

Full_FNs = cat(1,Full_FNs{:});

% Select the appropriate group
groupID = groupPFCCellsByDepth(Full_FNs,anatGroup);
SelectedFNs = Full_FNs(groupID);

[Location, Depth, Chan, Dist] = FindUnitLocation_mPFC(SelectedFNs);


if plotData
    nCellsTotal = size(Location,1);
    
    RatID = cellfun(@(x) x(1:4),SelectedFNs,'UniformOutput',0);
    AllRats = unique(RatID,'stable');
    colorOptions = jet(length(AllRats));
    
    fullColors = cellfun(@(x) colorOptions(ismember(AllRats,x),:)',RatID,'UniformOutput',0);
    fullColors = cat(2,fullColors{:})';
    
    RandAdd = (rand(nCellsTotal,1)-.5)/50;
    dummyIdx = cellfun(@(x) FastFind(ismember(RatID,x),1),AllRats);
    
   
    fh1 = figure;
    hold on 

    % Plot the cell locations
    for iD = 1:length(dummyIdx)
        idx = dummyIdx(iD);
        scatter3(Location(idx,1)+RandAdd(idx),Location(idx,2)+RandAdd(idx),-Location(idx,3)+RandAdd(idx),5,fullColors(idx,:));
    end
    Location(dummyIdx,:) = [];
    fullColors(dummyIdx,:) = [];
    RandAdd(dummyIdx) = [];
    
    scatter3(Location(:,1)+RandAdd,Location(:,2)+RandAdd,-Location(:,3)+RandAdd,5,fullColors);
    
     % Load & plot the anatomy subregion vertices
    PFCVertex = ExtractPFCAnatVrtx;
    fh1 = PlotAtlasPFCSubregions(fh1,PFCVertex,'plotType',anatplotType);
    figure(fh1)
    
    xlabel('ML')
    ylabel('AP')
    zlabel('DV')
    axis(axisVals)
%     plotXPlane(0)
%     legend(AllRats)
end
