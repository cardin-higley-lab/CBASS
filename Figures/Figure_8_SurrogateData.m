%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
iExp = 2;
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/'; chDir
% = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath, 'Figures', 'Figure_8_SurrogateData');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDir);
sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Test the method of surrogate data of different size.
% Set LFP and surrogate exemple 
inAnchor    = 40 * sREC.inSampleRate;
db1WinSec   = [0 6];

% Plots the figure
clear hFIG
hFIG = figure('Position', [50 50 1250 600]);
chFigName = chExperiment;

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

CBASS_SaveFig(chOutPath, hFIG, {chFigName});