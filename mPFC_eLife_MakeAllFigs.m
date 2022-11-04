function mPFC_eLife_MakeAllFigs

figOutputDir = 'E:\Data For Upload\mPFC Subregions - For Upload\Fig Outputs\';


%% Behavioral Data

mPFC_RevisionAnalysis_Behavior_eLife('figOutputDir',figOutputDir)


%% Recordings Data

mPFC_RevisionAnalysis_Recordings_eLife('figOutputDir',figOutputDir)


%% Example Cells

mPFC_RevisionAnalysis_Examples_eLife('figOutputDir',figOutputDir)


%% Regression of single cells to behavior

mPFC_RevisionAnalysis_SingleCellRegress_eLife('figOutputDir',figOutputDir)


%% Transfer Entropy to define Subregions

mPFC_RevisionAnalysis_TransferEntropy_eLife('figOutputDir',figOutputDir)


%% General Task Varaibles (Single Cell Regression)

mPFC_RevisionAnalysis_GeneralTaskVar_eLife('figOutputDir',figOutputDir)


%% Offer Zone Activitiy

mPFC_RevisionAnalysis_OfferZone_eLife('figOutputDir',figOutputDir)


%% Wait Zone Activity

mPFC_RevisionAnalysis_WaitZone_eLife('figOutputDir',figOutputDir)


%% Linger Zone Activity

mPFC_RevisionAnalysis_LingerZone_eLife('figOutputDir',figOutputDir)

