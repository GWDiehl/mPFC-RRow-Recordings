function [unitGroupings, stdLocation] = groupPFCCellsByDepth(S_FNs,dataGrouping,varargin)

% General function for grouping mPFC cells based on anatomical grouping
% seen by binning TC similarity by depth. Currently grouping is only done
% based on depth, or data that came from probes in VO. Depths are derived
% from rough quantification of similarity over anat depth plots.
%
% Inputs
% S_FNs - File names of the cells to be grouped.
% dataGrouping - Cell array of classes that you would like membership in.
%
% Ouptut unitGroupings - Cells x Groups logical identifying membership in
%       the requested group.

% GWD June 2020

VORats = {'R530_Si01' 'R535_Si01' 'R530_SiA1' 'R535_SiA1'};

PFCVertex = ExtractPFCAnatVrtx;
includeEdges = 1;

DVAnatBounds = [3.5 4.7];
DVAnatBounds_V2 = [3.8 4.7];
TE_AnatBounds = [2.9 3.85 4.45];

process_varargin(varargin);

stdLocation = FindUnitLocation_mPFC(S_FNs);
stdLocation = abs(stdLocation); % Make everything on one side

nCells = length(S_FNs);

% Find the data that came from the VO probes
unitProbe = cellfun(@(x) x(end-6:end-3),S_FNs,'UniformOutput',0);
unitRat = cellfun(@(x) x(1:4),S_FNs,'UniformOutput',0);
unitID = cellfun(@(x,y) [x,'_',y],unitRat,unitProbe,'UniformOutput',0);
VO_Cells = cellfun(@(x) ismember(x,VORats),unitID);

if ~iscell(dataGrouping)
    dataGrouping = {dataGrouping};
end

unitGroupings = false(nCells,length(dataGrouping));
for iG = 1:length(dataGrouping)
    
    % Start with everything as a valid candidate. Retain data according to
    % AND logic.
    validCells = true(nCells,1);
    
    excludeData = VO_Cells;
    
    targetAnat = cell(1);
    switch dataGrouping{iG}
        case {'All'}
            targetAnat{1} = [0 7];
            excludeData = false(size(excludeData));
        case {'None'}
            targetAnat{1} = [0 7];
            excludeData = true(size(excludeData));
            
            % Grouping derived based on DV location
        case {'GroupD' 'GroupA' 'Group1'}
            targetAnat{1} = [0 DVAnatBounds(1)];
        case {'GroupM' 'GroupB' 'Group2'}
            targetAnat{1} = DVAnatBounds;
        case {'GroupV' 'GroupC' 'Group3'}
            targetAnat{1} = [DVAnatBounds(2) 7];
        case {'GroupVO' 'Group4' 'VO'}
            targetAnat{1} = [0 7];
            excludeData = ~VO_Cells;
        case {'Group~VO' '~VO' 'No_VO'}
            targetAnat{1} = [0 7];
          
            % Grouping derived based on 3d location with respect to the
            % atlas
        case {'M2' 'M2_Atlas'}
            targetAnat{1} = PFCVertex.M2;
        case {'Acc' 'Acc_Atlas'}
            targetAnat{1} = PFCVertex.Acc;          
        case {'PL' 'PL_Atlas'}            
            targetAnat{1} = PFCVertex.PL;   
        case {'IL' 'IL_Atlas'}
            targetAnat{1} = PFCVertex.IL;  
        case {'MO' 'MO_Atlas'}
            targetAnat{1} = PFCVertex.MO;   
        case {'DP' 'DP_Atlas'}
            targetAnat{1} = PFCVertex.DP;   
            
            % Grouping derived as a combination of factors
        case {'dPL' 'PL_Group1' 'PL_Dorsal'} % Dorsal half of PL
            targetAnat = cell(1,2);
            targetAnat{1} = PFCVertex.PL;
            targetAnat{2} = [0 DVAnatBounds(1)]; 
            
        case {'vPL' 'PL_Group2' 'PL_Ventral'}  % Ventral half of PL
            targetAnat = cell(1,2);          
            targetAnat{1} = PFCVertex.PL;
            targetAnat{2} = [DVAnatBounds(1) 7];         
            
            % Use the other dPL/vPL boundary at about 3.8mm
        case {'dPL_V2'} % Dorsal half of PL
            targetAnat = cell(1,2);
            targetAnat{1} = PFCVertex.PL;
            targetAnat{2} = [0 DVAnatBounds_V2(1)]; 
            
        case {'vPL_V2'} % Ventral half of PL
            targetAnat = cell(1,2);
            targetAnat{1} = PFCVertex.PL;
            targetAnat{2} = [DVAnatBounds_V2(1) 7];
            
            % Four mPFC subregions based on the boundaries extracted from
            % TE analysis
        case {'Acc_TE'}
            targetAnat{1} = [0 TE_AnatBounds(1)];
        case {'dPL_TE'}
            targetAnat{1} = [TE_AnatBounds(1) TE_AnatBounds(2)];
        case {'vPL_TE'}
            targetAnat{1} = [TE_AnatBounds(2) TE_AnatBounds(3)];
        case {'IL_TE'}
            targetAnat{1} = [TE_AnatBounds(3) 7];
            
        case {'dmPFC' 'dorsalPFC'}
            targetAnat{1} = [0 TE_AnatBounds(2)];
        case {'vmPFC' 'ventralPFC'}
            targetAnat{1} = [TE_AnatBounds(2) 7];
            
            
        otherwise % Treat as none
            targetAnat{1} = [0 7];
            excludeData = true(size(excludeData));
    end
    % Remove any raw exclusions
    validCells(excludeData) = 0;
    
    % For each targetAnat entry, loop through and keep only the applicable
    % cells that match each and every criteria (AND logic)
    for iT = 1:length(targetAnat)    
        if numel(targetAnat{iT}) == 2
            % Find cells @ depth
            currCells = stdLocation(:,3) >= targetAnat{iT}(1) & stdLocation(:,3) <= targetAnat{iT}(2);
            
        else
            % Identify based on atlas derived vertices
            
            % Grab all of the available AP section levels
            APLevels = unique(PFCVertex.External(:,2));
            
            currCells = false(nCells,1);
            for iC = 1:nCells
                
                % Find the section that most closly match the cell location
                [~,minIdx] = min(abs(stdLocation(iC,2) - APLevels));
                validTargetAnat = targetAnat{iT}(:,2) == APLevels(minIdx);
                % If the targetAnat group does not have the required section it
                % does not qualify, move along.
                if ~any(validTargetAnat)
                    continue
                end
                
                % Is the cell location within the polygon of the closest
                % section
                [inTA, onTA] = inpolygon(stdLocation(iC,1),stdLocation(iC,3),targetAnat{iT}(validTargetAnat(:),1),targetAnat{iT}(validTargetAnat(:),3));
                % Does on the edge count?
                if includeEdges
                    currCells(iC) = inTA||onTA;
                else
                    currCells(iC) = inTA;
                end
            end
            
        end
        % Valid cells are those that have continued to meet all criteria
        validCells = validCells & currCells;
    end
    
    % Logical those cells that are valid
    unitGroupings(validCells,iG) = 1;
    
end