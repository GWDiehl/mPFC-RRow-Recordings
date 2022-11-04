function mPFC_eLife_RunAllAnalyses

%% Averge FR in each maze zone per lap

mPFC_ExtractOverallZoneRates

%% Linear Regression between cell FR and task variables

mPFC_LinearRegressionToBehavior

%% Transfer Entropy between Cells

mPFC_ComputeByChunkTE

%% FR PETHs for each maze zone

mPFC_PETHsAcrossZones


%% MI for offer Zone 

mPFC_OfferZoneMI

%% MI for Wait zone

mPFC_WaitZoneMI

