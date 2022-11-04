function  S = loadAllSpikes_GWD(S_FNs,varargin)

restrictToTaskTime = 0;
TFileRoot = '.t';

[RootDir, ~, ~, ~, ~, ~, SpikesDir] = ID_GWD_Dirs;


process_varargin(varargin);


RatID = cellfun(@(x) x(1:4),S_FNs,'UniformOutput',0);
SSN = cellfun(@(x) x(1:15),S_FNs,'UniformOutput',0);


FullSpikeFile = cellfun(@(x,y,z) fullfile(RootDir,SpikesDir,y,z,[x TFileRoot]),S_FNs,RatID,SSN,'UniformOutput',0);
S = LoadSpikes(FullSpikeFile);


if restrictToTaskTime
    origDir = pwd;
    uniqueSSNs = unique(SSN);
    
    for iS = 1:length(uniqueSSNs)    
        fullIdx = ismember(SSN,uniqueSSNs{iS});
        
        PathFN = fullfile(RootDir,BehavDir,uniqueSSNs{iS}(1:4),uniqueSSNs{iS},[uniqueSSNs{iS} PathFileRoot]);
        keysFN = [strrep(uniqueSSNs{iS},'-','_') KeysFileRoot];       
        
        cd(fileparts(PathFN))
        [~,keysFunc] = fileparts(keysFN);
        keys = evalKeysFile(keysFunc);
        
        startTask = keys.TimeOnTrack;
        endTask = keys.TimeOffTrack; 
        
        S(fullIdx) = cellfun(@(x) x.restrict(startTask,endTask),S(fullIdx),'uniformoutput',0);
    end
    
    cd(origDir);
end




