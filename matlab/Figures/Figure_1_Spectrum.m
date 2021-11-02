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
chOutPath   = fullfile(chOutDir, 'Figure_1_Spectrum');

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDataPath);
% sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Sets general parameter
% Sets the band and the state of interest
chLabel         = 'Gamma';
db1Band         = [30 80];
bl1State        = sREC.bl1Run;
chStateLabel    = 'Running';

% Sets options for trough extraction
chDataFormat    = 'complex';
blZScore        = true;
inRefChan       = 5;

% Sets parameters for plotting the trough
in1StateON  = find(~bl1State(1:end - 1) & bl1State(2:end));
inAnchor    = in1StateON(2) - (0.2 * sREC.inSampleRate);
db1WinSec   = [0 .6];

% Sets the gloabal label for the experiment;
chExpLbl    = [cEXP{iExp} '_' chLabel];

% Gets the troughs
sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
%% Initializizes the figure
clear hFIG
hFIG = figure('Position', [50 50 1250 600]);
cFIG_NAME = {chExpLbl};

% Plots the spectrum
CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate, false, bl1State); hold on
legend('AutoUpdate', 'off');
db1YL = ylim; plot(db1Band(1) * [1 1], db1YL, '--r'); plot(db1Band(2) * [1 1], db1YL, '--r');
xlim([0 120])
title([chLabel ' Power vs ' chStateLabel]);

%% Saves the figure
CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
close all;