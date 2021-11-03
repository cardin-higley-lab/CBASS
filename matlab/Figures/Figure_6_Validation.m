%% Loads data
clear, close all

%Adds CBASS to the path
chCBASSDir  = 'D:\CBASS\matlab\';
addpath(genpath(chCBASSDir)); %Add the root folder to current folder

%Sets path to the data
% chDataDir   = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/data';
chDataDir   = 'D:\gamma_bouts\';
cEXP        = {'Example_1', 'Example_2'};
iExp        = 1;
chDataPath  = fullfile(chDataDir, cEXP{iExp});

%Sets the output path
chOutDir    = 'D:\CBASS\matlab\Figures';
chOutPath   = fullfile(chOutDir, 'Figure_6_Validation');

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDataPath);
%% Sets options
cBAND           = {[30 80]};
cSTATE          = {sREC.bl1Run};

sOPTION.cBAND_LABEL     = {'Gamma'};  % character array or cell array of labels  for bands in cBAND (i.e. {'Beta', 'Gamma'})
sOPTION.cSTATE_LBL      = {'Running'};  % character array or cell array of labels  for states in cSTATE (i.e. {'Stim', 'Running'})
sOPTION.blVerbose       = true;    % Sets whether to print progress on the screen. Default is true 

% L1 options for formatting hilbert troughs (Function CBASS_L1_GetTrough)
sOPTION.chDataFormat    = 'complex'; % Format of hilbert transform coordinates ('polar' or 'complex' Default is 'complex')
sOPTION.inRefChan       = 5;    % Reference channel (Default is the last of row of db2LFP.

% Runs CBASS
sFREQ_BAND = CBASS_Main_DetectEvents(sREC.db2LFP, sREC.inSampleRate, cBAND, cSTATE, sOPTION);

%% Plots the figure
clear hFIG
hFIG = figure('Position', [50 50 1250 600]);
cFIG_NAME = {[cEXP{iExp} '_' sOPTION.cSTATE_LBL{1}]};

% Plots CSD 
subplot(1, 3, 1)
CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sFREQ_BAND.bl1Event, [-1 1] ./ mean(cBAND{1}), sREC.inSampleRate);
title(sprintf('CSD around Events'));

% Plots event enrichment 
subplot(1, 3, 2)
chSigString = CBASS_Plot_Pulse_vs_State(sFREQ_BAND.bl1Event, cSTATE{1}, [-5 10], sREC.inSampleRate);
title(sprintf('Event Rate around %s Onset', sOPTION.cSTATE_LBL{1}))

% Plots the power under and above the 80%percentile of a moving average of
% the number of pulses
dbWinLenSec     = 0.5;
db1Win          = rectwin(dbWinLenSec * sREC.inSampleRate);
db1NPulse       = conv(sFREQ_BAND.bl1Event, db1Win, 'same')./dbWinLenSec;
db80pct         = quantile(db1NPulse, .8);
bl180pctPlus    = db1NPulse > db80pct;
subplot(1, 3, 3);
CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate,...
    false, bl180pctPlus, {sprintf('<%.1fHz', db80pct) sprintf('>%.1fHz', db80pct)});
title('LFP Power vs Event Frequency')

CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME, 'png');
close all;