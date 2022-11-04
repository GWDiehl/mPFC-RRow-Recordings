function mPFC_RevisionAnalysis_Recordings_eLife(varargin)

figOutputDir = 'E:\DATA-Geoff\GDiehl Results\mPFC_RRow Manuscript\eLife\';

process_varargin(varargin);


CleanMatlabPlotDefaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   General Recordings   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


currSubDir = 'Recordings\';

if ~exist([figOutputDir,currSubDir],'dir')
    mkdir([figOutputDir,currSubDir]);
end

%% Example WF

[RootDir, ~, ~, ~, ~, ~, ~, ~, ChannelMapDir] ...
    = ID_GWD_Dirs;

S_FN = 'R535-2019-07-27-Si02_06';
chanMap = [RootDir,ChannelMapDir,'\R535_channel_map_Si_ProbeC.txt'];

setJetAsDefaultColorMap
[fh1, fh2,sortedWF,peakChan] = PlotSortedWF(S_FN,chanMap);

figFN = [figOutputDir,currSubDir,'ExampleWF'];
saveas(fh1,figFN,'svg'); close(fh1)

figFN = [figOutputDir,currSubDir,'ExampleWF_Image'];
saveas(fh2,figFN,'svg'); close(fh2)


%% 3D localization of cells

[~, ~, ~, ~, ~,fh] = LocalizeUnits_mPFC('anatGroup','No_VO','anatplotType','Both');

figure(fh)
view(73,36)
figFN = [figOutputDir,currSubDir,'All3DLoc'];
saveas(fh,figFN,'fig');
saveas(fh,figFN,'svg');

view(0,12)
figFN = [figOutputDir,currSubDir,'All3DLocV2'];
saveas(fh,figFN,'svg');

view(91,4)
figFN = [figOutputDir,currSubDir,'All3DLocV3'];
saveas(fh,figFN,'svg');

close(fh)

% Make a rotating GIF

movieOut = [figOutputDir,currSubDir];
movieFN = '3DAnatWCells.gif';
fh = MakePFCAnatMovie(movieOut,movieFN);
close(fh)

%% Breakdown of cell recordings per rat/subregion etc.

% Prelim Directoires & DataDef
[rootDir, DataDefName] = ID_GWD_Dirs;
dataDef = load([rootDir,DataDefName]);
OI = ismember(dataDef.Task,'RRow') & ismember(dataDef.ExpType,'Recording-Stable');

[SSNs,~,~,~,S_FNs] = GWD_LoadAllRRowSessions('dataToLoad',{'Spikes'});
RatID = cellfun(@(x) x(1:4),SSNs,'UniformOutput',0);
allRats = unique(RatID,'stable');

% Remove any VO recordings from the cell list
no_VO = cellfun(@(x) groupPFCCellsByDepth(x,{'No_VO'}),S_FNs,'UniformOutput',0);
S_FNs = cellfun(@(x,y) x(y),S_FNs,no_VO,'UniformOutput',0);
cellCounts = cellfun(@length,S_FNs);

tableFields = {'RatID' 'TotalSess' 'TotalCells' 'AvgCells' 'CellsRange' 'CellsStd' 'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE' 'Acc_TEAvg' 'dPL_TEAvg' 'vPL_TEAvg' 'IL_TEAvg'};
dataCountsfigOutputDir = struct;
for iT = 1:length(tableFields)
    switch tableFields{iT}
        case 'RatID'
            currData = allRats;
        case 'TotalSess'
            currData = cellfun(@(x) sum(ismember(dataDef.RatID(OI),x)),allRats);
        case 'TotalCells'
            currData = cellfun(@(x) sum(cellCounts(ismember(RatID,x))),allRats);
        case 'AvgCells'
            currData = cellfun(@(x) nanmean(cellCounts(ismember(RatID,x))),allRats);
        case 'CellsRange'
            currData = cellfun(@(x) min(cellCounts(ismember(RatID,x))),allRats);
            currData(:,2) = cellfun(@(x) max(cellCounts(ismember(RatID,x))),allRats);
        case 'CellsStd'
            currData = cellfun(@(x) nanstd(cellCounts(ismember(RatID,x))),allRats);
            
        case {'Acc' 'dPL_V2' 'vPL_V2' 'IL' 'AccAvg' 'dPL_V2Avg' 'vPL_V2Avg' 'ILAvg' 'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE' 'Acc_TEAvg' 'dPL_TEAvg' 'vPL_TEAvg' 'IL_TEAvg'}
            cellGroup = cellfun(@(x) groupPFCCellsByDepth(x,{strrep(tableFields{iT},'Avg','')}),S_FNs,'UniformOutput',0);
            grp_S_FNs = cellfun(@(x,y) x(y),S_FNs,cellGroup,'UniformOutput',0);
            grpCount = cellfun(@length,grp_S_FNs);
            if contains(tableFields{iT},'Avg')
                currData = cellfun(@(x) nanmean(grpCount(ismember(RatID,x))),allRats);
            else
                currData = cellfun(@(x) sum(grpCount(ismember(RatID,x))),allRats);
            end
    end
    dataCountsfigOutputDir.(tableFields{iT}) = currData';
end
outFN = [figOutputDir,currSubDir,'CountsfigOutputDir'];
save(outFN,'dataCountsfigOutputDir');



%% Identify Principal cells vs Interneurons

% Done based on grouping data according to avg firing rate and spike width.
% Addition of PeakValleyRatio as a criteria does not emperically add
% anything in this case. Clustering via KMeans or a Gausian mixed model
% both dont do well so draw the boundary by hand and log the results.

[RootDir, ~, ~, ~, ~, ~, ~, ~, ~, AnalysisDir] = ID_GWD_Dirs;
outputFN = fullfile(RootDir,AnalysisDir,'CellTypeClassification');
spikingQuants = {'AvgRate','SpikeWidth','PeakValleyRatio'};
anatGroupings = {'No_VO' 'Acc_TE' 'dPL_TE' 'vPL_TE' 'IL_TE'};
saveNewResults = 0;
quantBounds = {-7.5:2.5:72.5,0:.05:.9};


[~,~,~,S,S_FNs] = GWD_LoadAllRRowSessions('dataToLoad',{'Spikes'});
S = cat(1,S{:});
S_FNs = cat(1,S_FNs{:});
validCells = groupPFCCellsByDepth(S_FNs,'No_VO');

[fhScatter,fhBox,spikeQuantData] = BasicSpikingVsAnat_V2(S,S_FNs,spikingQuants);

for iQ = 1:length(spikingQuants)
    close(fhScatter{iQ})
    close(fhBox{iQ})
end

% Plot the cell parameters and draw a boundary between Principal and
% Interneurons
if saveNewResults
    fh1 = figure;
    scatter(spikeQuantData(:,1),spikeQuantData(:,1))
    xlabel(spikingQuants{1})
    ylabel(spikingQuants{2})
    [x,y] = ginput;
    close(fh1)
else
    pastGrouping = load(outputFN);
    x = pastGrouping.BoundaryPoints(:,1);
    y = pastGrouping.BoundaryPoints(:,2);
end

idx = double(inpolygon(spikeQuantData(:,1)',spikeQuantData(:,2)',x,y));
idx(~idx) = 2;

fh1 = figure; 
subplot(1,2,1); hold on
scatter(spikeQuantData(idx == 1,1),spikeQuantData(idx == 1,2))
scatter(spikeQuantData(idx == 2,1),spikeQuantData(idx == 2,2))
plot(x,y)
xlabel(spikingQuants{1})
ylabel(spikingQuants{2})
xlim([quantBounds{1}(1) quantBounds{1}(end)])
ylim([quantBounds{2}(1) quantBounds{2}(end)])
legend({'Principal' 'Interneuron' 'Boundary'})

subplot(1,2,2); hold on
cellCount = histcn(spikeQuantData(:,1:2), quantBounds{:});
s = imagesc(computeBinCenters(quantBounds{1}),computeBinCenters(quantBounds{2}),log10(cellCount)');
set(s,'alphadata',cellCount' >0)
box off
xlim([quantBounds{1}(1) quantBounds{1}(end)])
ylim([quantBounds{2}(1) quantBounds{2}(end)])
axis xy
xlabel(spikingQuants{1})
ylabel(spikingQuants{2})
caxis([0 2.5])
plot(x,y,'r')
colorbar

% Save the classification results and identified boundary
CellClass = struct;
CellClass.CellID = S_FNs;
CellClass.Interneuron = idx == 1;
CellClass.Principal = idx == 2;
CellClass.BoundaryPoints = cat(2,x,y);

if saveNewResults
    save(outputFN,'-struct','CellClass')
    
    anatGrpFN = fullfile(RootDir,AnalysisDir,'CellAnatGroupings');
    ExtractAndSavePFCCellLoc(S_FNs,anatGroupings,anatGrpFN)
end

figFN = [figOutputDir,currSubDir,'Principal-InterneuronSplit'];
fh1.Position = [100 100 1150 620];
saveas(fh1,figFN,'svg'); close(fh1)


%% General firing properties as a function of depth

spikingQuants = {'AvgRate' 'MedianISI' 'ISIStd' 'ISIRatio' 'COV' 'SpikeWidth' 'PeakValleyRatio' 'ISIRatio_Log'}; %

cellGroupings = {'Principal' 'Interneuron' 'No_VO'};
groupNames = {'Principal' 'Interneuron' 'All'};
plotScatter = [1 1 0];
plotErrorbar = [0 0 1];

[~,~,~,S,S_FNs] = GWD_LoadAllRRowSessions('dataToLoad',{'Spikes'});
S = cat(1,S{:});
S_FNs = cat(1,S_FNs{:});

[validCells,stdLoc] = groupPFCCellsByDepth(S_FNs,'No_VO');
stdDepth = stdLoc(validCells,3);

[fhScatter,fhBox,fullData] = BasicSpikingVsAnat_V2(S(validCells),S_FNs(validCells),spikingQuants,...
    'cellGroupings',cellGroupings,'groupNames',groupNames,'plotScatter',plotScatter,'plotErrorbar',plotErrorbar);

figure(fhScatter{ismember(spikingQuants,'ISIRatio')})
set(gca, 'YScale', 'log')
ylim([0 200])


corr_Coeff = nan(length(spikingQuants),length(cellGroupings));
corr_pVal = nan(length(spikingQuants),length(cellGroupings));
slope = nan(length(spikingQuants),length(cellGroupings));

for iQ = 1:length(spikingQuants)
    figFN = [figOutputDir,currSubDir,spikingQuants{iQ}];
    saveas(fhScatter{iQ},figFN,'svg'); close(fhScatter{iQ})
    
    %     figFN = [figOutDir,currSubDir,spikingQuants{iQ},'_Boxplot'];
    %     saveas(fhBox{iQ},figFN,'svg');
    close(fhBox{iQ})
    
    for iG = 1:length(cellGroupings)
        valid = isfinite(fullData(:,iQ,iG));
        [corr_Coeff(iQ,iG), corr_pVal(iQ,iG)] = nancorr(stdDepth(valid),fullData(valid,iQ,iG));
        temp = polyfit(stdDepth(valid),fullData(valid,iQ,iG),1);
        slope(iQ,iG) = temp(1);
    end
end





