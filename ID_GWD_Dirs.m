function [RootDir, DataDefName, BehaviorDir, TaskSuffix, PathSuffix, KeysSuffix, UnitsDir, LFPDir, ChannelMapDir, AnalysisDir, PromotionDir] ...
    = ID_GWD_Dirs(varargin) 
    
% Load up the standard directories of where data live and are names for GWD

% For use on other systems, by other users, or with data that has been
% moved/copied, etc. these fields should be changed according to the
% location of the stored data. Within each of these parents, data/analyses
% are in standard byRat\bySession\DATA organization


% This root Dir will NEED to be updated according to your own needs
RootDir = 'E:\Data For Upload\mPFC Subregions - For Upload\Data\';

% Everything else, to each their own
BehaviorDir = 'Processed Behavior\';
AnalysisDir = 'Analysis\';
UnitsDir = 'Processed Phys\Processed Units\';
LFPDir = 'Processed Phys\Processed LFP\';
ChannelMapDir = 'Misc\Channel Maps\';

PromotionDir = 'Promotable Data\';

DataDefName = 'Misc\DataDef_RRow';

TaskSuffix = '-RRow.mat';
PathSuffix = '-vt.mat';
KeysSuffix = '_keys.m';

process_varargin(varargin);



