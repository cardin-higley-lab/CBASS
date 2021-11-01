%% Loads data
%clear all
%close all
chExperiment = 'Example_1';
chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
chDir       = fullfile(chPath,'data', chExperiment);
% chPath      = 'D:\gamma_bouts\';
% chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Matlab_Pipeline', 'output');
chTmpPath   = fullfile(chPath,'Matlab_Pipeline', 'tmp/');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir')
   mkdir(chTmpPath)
end
if ~exist(chOutPath, 'dir')
   mkdir(chOutPath)
end
%%
sREC        = CBASS_L0_LoadData(chDir);
sREC        = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Computes the troughs
chLabel         = 'Beta';
dbBand          = [15 30];
inRefChan       = 5;
chDataFormat    = 'complex';

sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
%% Call python to do the embeding  
chMethod        = 'umap';
blZScore        = true;
inN_Component   = 3;
inLabelTags     = ['any_labels'];
chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat];

[status, commandOut] = CBASS_L2_Embed(sTROUGH, chTmpPath, chOutPath, chExpLabel, ...
    chMethod, chDataFormat, blZScore, inN_Component)
%% Identifies enriched regions
[sREGION, sCLU] = CBASS_L2_GetEnrichedRegion(sTROUGH, sREC.bl1Pres, blZScore);
%% Generates the embeding
in1EmbedLabel    = sREGION(1).bl1Member';

[status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, inLabelTags, ...
    chMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel)
%% Get the pulses
sPULSE      = CBASS_L3_GetPulse(sREC, sTROUGH, sREGION(1).bl1Member);
%% Plots the histogram for the pulses
hFIG = figure;
chFigName = ['Peak_Histogram_' chExpLabel];
CBASS_Plot_PeakHistogram(sPULSE);
CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
%% Sets the event of interst and plots event triggered average
cEVENT          = {sPULSE.bl1Pulse sTROUGH.in1Index sTROUGH.in1Index(sREGION(1).bl1Member)};
cEVENT_LABEL    = {'Pulse', 'Everything', 'SelTrough'};
inNEvt          = length(cEVENT);

hFIG = figure('Position', [50 50 1250 600]);
chFigName = ['Triggered_Averages_' chExpLabel];
for iEvt = 1:inNEvt
    subplot(2, inNEvt, iEvt)
    CBASS_Plot_LFP_EventTriggeredAverage(sREC.db2LFP, cEVENT{iEvt}, [-.075 .075], sREC.inSampleRate);
    title(sprintf('%s Triggered Average: LFP', cEVENT_LABEL{iEvt}))
    subplot(2, inNEvt, inNEvt + iEvt)
    CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, cEVENT{iEvt}, [-.075 .075], sREC.inSampleRate);
    title('CSD')
end
CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
%% Plots State enrichment of pulses
cEVENT          = {sPULSE.bl1Pulse sTROUGH.in1Index(sREGION(1).bl1Member)};
cEVENT_LABEL    = {'Pulse', 'SelTrough'};
inNEvt          = length(cEVENT);

hFIG = figure('Position', [50 50 1250 600]);
chFigName = ['PulseVsState_' chExpLabel];
for iEvt = 1:inNEvt
    subplot(1, 2, iEvt)
    chSigString = CBASS_Plot_Pulse_vs_State(cEVENT{iEvt}, sREC.bl1Pres, [-5 10], sREC.inSampleRate);
    title(sprintf('%s : %s', cEVENT_LABEL{iEvt}, chSigString))
end
CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
%% Plots LFP example
in1Anchor       = find(~sREC.bl1Pres(1:end - 1) & sREC.bl1Pres(2:end));
cEVENT          = {sPULSE.bl1Pulse  sTROUGH.in1Index(sREGION(1).bl1Member)};
cEVENT_LABEL    = {'Pulse', 'SelTrough'};

hFIG = figure('Position', [50 50 1250 600]);
chFigName = ['LFP_Examples_' chExpLabel];
for iEx = 1:4
    subplot(2, 2, iEx)
    CBASS_Plot_LFP_EventExample(sREC.db2LFP, in1Anchor(iEx), [-2 3], sREC.inSampleRate, cEVENT, cEVENT_LABEL)
    title(sprintf('Data around %s onset: %d', 'Stim', iEx));
end
CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');