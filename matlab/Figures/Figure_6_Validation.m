%% Loads data
clear, close all

%Adds CBASS to the path
chCBASSDir  = 'D:\CBASS\matlab\';
addpath(genpath(chCBASSDir)); %Add the root folder to current folder

%Sets path to the data
% chDataDir   = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/data';
chDataDir   = 'D:\gamma_bouts\';
cEXP        = {'Example_1', 'Example_2'};
iExp        = 2;
chDataPath  = fullfile(chDataDir, cEXP{iExp});

%Sets the output path
chOutDir    = 'D:\CBASS\matlab\Figures';
chOutPath   = fullfile(chOutDir, 'Figure_9_Validation');

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDataPath);
sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);

% Sets general parameter
chLabel         = 'Beta';
db1Band         = [15 30];
bl1State        = sREC.bl1Pres;
chStateLabel    = 'Stim';

% Sets options for enriched region identification
blUseRate   = false;
inMethod    = 1;
inNClu      = 20;

% Sets options for trough extraction
chDataFormat    = 'complex';
blZScore        = true;
inRefChan       = 5;

% Sets the gloabal label for the experiment;
chExpLbl    = [chExperiment '_' chLabel];

% Gets the troughs
sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Gets the rand trough
sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Identifies enriched regions
sFILTER     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, bl1State, blZScore, 1, inNClu);

% Get the pulses for each enriched region
sPULSE      = CBASS_L3_GetPulse_2(sREC, sFILTER(1).db2Filter, true, true);

% Initializizes the figure
clear hFIG
hFIG = figure('Position', [50 50 1250 600]);
cFIG_NAME = {chExpLbl};

% Plots CSD 
subplot(1, 3, 1)
CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sPULSE(1).bl1Pulse, [-1 1] ./ mean(db1Band), sREC.inSampleRate);
title(sprintf('CSD around pulse'));

% Plots event enrichment 
subplot(1, 3, 2)
chSigString = CBASS_Plot_Pulse_vs_State(sPULSE.bl1Pulse, bl1State, [-5 10], sREC.inSampleRate);
title(sprintf('%s', chSigString))

% Plots the power under and above the 80%percentile of a moving average of
% the number of pulses
dbWinLenSec     = 0.5;
db1Win          = rectwin(dbWinLenSec * sREC.inSampleRate);
db1NPulse       = conv(sPULSE.bl1Pulse, db1Win, 'same')./dbWinLenSec;
db80pct         = quantile(db1NPulse, .8);
bl180pctPlus    = db1NPulse > db80pct;
subplot(1, 3, 3);
CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate,...
    false, bl180pctPlus, {sprintf('<%.1fHz', db80pct) sprintf('>%.1fHz', db80pct)});
title([chLabel ' Power vs Pulse Frequency'])

CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
close all;