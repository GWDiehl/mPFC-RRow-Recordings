
function fh = PlotAtlasPFCSubregions(fh,pFCVertex,varargin)


% Takes a structure identifying the vertices of subregion bounds defined by
% the paxinos & watson atlas and plots them on top of an existing figure.
% If no figure is input ([]) it will create a new one.
%
% Specific subregions can be excluded using the subregionsToExlude flag.

% GWD March 2021


subregionsToExclude = {'External'};
colorScheme = 'Jet';

plotType = 'Sections';

opacityRange = [.05 .15];

process_varargin(varargin);


%%

if isempty(fh)
    fh = figure;
end
figure(fh)
hold on

% Find all of the AP Levels in the data and which subregions are listed
APLevels = unique(pFCVertex.External(:,2));
subregions = fieldnames(pFCVertex);

% Remove any unwanted subregions
for iE = 1:length(subregionsToExclude)
    subregions(ismember(subregions,subregionsToExclude{iE})) = [];
end
nSubregions = length(subregions);

opacityVals = linspace(opacityRange(1),opacityRange(2),length(APLevels)-1);


% Identify the appropriate color scheme. Default Jet
if ischar(colorScheme)
    switch colorScheme
        case {'Jet' 'jet'}
            sectionColors = jet(nSubregions);
        otherwise
            error('Color scheme not yet implemented')
    end
elseif iscell(colorScheme)
    sectionColors = cat(1,colorScheme{:});
end
assert(size(sectionColors,2)==3 & size(sectionColors,1)>=nSubregions,'Your color scheme will not work for plotting')


% Plot the anatomy
for iE = 1:nSubregions
    currVtx = pFCVertex.(subregions{iE});
    currVtx(:,3) = -currVtx(:,3);
    
    switch plotType
        case 'Sections'
            for iS = 1:length(APLevels)
                validSlice = currVtx(:,2) == APLevels(iS);
                if ~any(validSlice)
                    continue
                end
                plot3(currVtx(validSlice,1),repmat(APLevels(iS),sum(validSlice),1),currVtx(validSlice,3),'color',sectionColors(iE,:));
                plot3(-currVtx(validSlice,1),repmat(APLevels(iS),sum(validSlice),1),currVtx(validSlice,3),'color',sectionColors(iE,:));
            end
            
        case 'Filled'
            for iS = 1:length(APLevels)-1
                validSlice = currVtx(:,2) == APLevels(iS);
                nextSlice = currVtx(:,2) == APLevels(iS+1);
                if ~any(validSlice) || ~any(nextSlice)
                    continue
                end
                
                biSecVertex = unique(cat(1,currVtx(validSlice,:),currVtx(nextSlice,:)),'rows');
                DT = delaunayTriangulation(biSecVertex(:,1),biSecVertex(:,2),biSecVertex(:,3));
                [F,P] = freeBoundary(DT);
                
                trisurf(F,P(:,1),P(:,2),P(:,3),...
                    'facecolor',sectionColors(iE,:),'FaceAlpha',opacityVals(iS),'edgecolor','none');
                trisurf(F,-P(:,1),P(:,2),P(:,3),...
                    'facecolor',sectionColors(iE,:),'FaceAlpha',opacityVals(end),'edgecolor','none');
            end
            
        case 'Both'
            for iS = 1:length(APLevels)
                validSlice = currVtx(:,2) == APLevels(iS);
                if ~any(validSlice)
                    continue
                end
                plot3(currVtx(validSlice,1),repmat(APLevels(iS),sum(validSlice),1),currVtx(validSlice,3),'color',sectionColors(iE,:));
                plot3(-currVtx(validSlice,1),repmat(APLevels(iS),sum(validSlice),1),currVtx(validSlice,3),'color',sectionColors(iE,:));
                
                if iS == length(APLevels)
                    continue
                end
                nextSlice = currVtx(:,2) == APLevels(iS+1);
                if ~any(nextSlice)
                    continue
                end
                
                biSecVertex = unique(cat(1,currVtx(validSlice,:),currVtx(nextSlice,:)),'rows');
                DT = delaunayTriangulation(biSecVertex(:,1),biSecVertex(:,2),biSecVertex(:,3));
                [F,P] = freeBoundary(DT);
                
                trisurf(F,P(:,1),P(:,2),P(:,3),...
                    'facecolor',sectionColors(iE,:),'FaceAlpha',opacityVals(iS),'edgecolor','none');
                trisurf(F,-P(:,1),P(:,2),P(:,3),...
                    'facecolor',sectionColors(iE,:),'FaceAlpha',opacityVals(end),'edgecolor','none');
                
            end
    end
end


