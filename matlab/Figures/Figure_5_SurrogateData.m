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
chOutPath   = fullfile(chOutDir,  'Figure_5_SurrogateData');

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDataPath);
sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Test the method of surrogate data of different size.
% Set LFP and surrogate exemple 
inAnchor    = 50 * sREC.inSampleRate;
db1WinSec   = [0 6];

% Plots the figure
clear hFIG
hFIG = figure('Position', [50 50 900 500]);
chFigName =  cEXP{iExp};

% Plots the exemples
subplot(2, 2, 1); 
CBASS_Plot_LFP_EventExample(sREC.db2LFP, inAnchor, db1WinSec, sREC.inSampleRate);
title('LFP');
subplot(2, 2, 2); 
CBASS_Plot_LFP_EventExample(sREC.db2LFP_Rnd, inAnchor, db1WinSec, sREC.inSampleRate);
title(sprintf('Surrogate LFP'));

% Plots the LFP extracts for real and surrogate data
clear hPLT
hPLT(1) = subplot(2, 2, 3);  CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate, true)
hPLT(2) = subplot(2, 2, 4);  CBASS_Plot_LFP_FourierPower(sREC.db2LFP_Rnd, sREC.inSampleRate, true)
linkaxes(hPLT, 'y');

CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
close all