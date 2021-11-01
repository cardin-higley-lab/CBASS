%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
% for iExp = 1:length(cEXP)
iExp = 1;
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Tests', 'Test_PowerSpectrum');
chTmpPath   = fullfile(chPath,'Tests', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDir);
sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
% Sets general parameter
%% Sets the loop
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND       = {[15 30], [30 80]};
cSTATE      = {sREC.bl1Pres sREC.bl1Run};
cSTATE_LBL  = {'Stim', 'Running'};

% Downsampels the LFP four time for the morlet transform
inFactor    = 4;
db2LFP_Tst  = [];
for iChan = 1:size(sREC.db2LFP, 1)
    db2LFP_Tst = cat(1, db2LFP_Tst, CBASS_U_DownSampleTrace(sREC.db2LFP(iChan, :), inFactor));
end

for iBnd = 1:length(cBAND)
    % Sets the gloabal label for the experiment;
    chLabel     = cSTATE_LBL{iBnd};
    chExpLbl    = [chExperiment '_' chLabel];
    
    % Initializizes the figure
    clear hFIG
    hFIG = figure('Position', [50 50 1250 600]);
    cFIG_NAME = {chExpLbl};
    
    % Dowsamples the the state vector for Wavelet spectra
    bl1State_DS = CBASS_U_DownSampleTrace(cSTATE{iBnd}, inFactor, 1);
    
    % Plots the Fourier spectra
        % Raw
    subplot(2, 3, 1);
    CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate, false, cSTATE{iBnd});
    title([cBAND_LABEL{iBnd} 'Welsh Power vs ' cSTATE_LBL{iBnd}]);
        
        % ZScored
    subplot(2, 3, 2);
    CBASS_Plot_LFP_FourierPower(zscore(sREC.db2LFP), sREC.inSampleRate, false, cSTATE{iBnd});
    title([cBAND_LABEL{iBnd} 'Welsh Zscore Power vs ' cSTATE_LBL{iBnd}]);
        
        % Bipolar Derivation
    subplot(2, 3, 3);
    CBASS_Plot_LFP_FourierPower(diff(sREC.db2LFP), sREC.inSampleRate, false, cSTATE{iBnd});
    title([cBAND_LABEL{iBnd} 'Welsh BiPol Power vs ' cSTATE_LBL{iBnd}]);
    
    % Plots the Morlet spectra
        % Raw
    subplot(2, 3, 4);
    CBASS_Plot_LFP_MorletPower(db2LFP_Tst, sREC.inSampleRate/inFactor, false, bl1State_DS);
    title([cBAND_LABEL{iBnd} 'Morlet Power vs ' cSTATE_LBL{iBnd}]);   
    
        % ZScored
    subplot(2, 3, 5);
    CBASS_Plot_LFP_MorletPower(zscore(db2LFP_Tst), sREC.inSampleRate/inFactor, false, bl1State_DS);
    title([cBAND_LABEL{iBnd} 'Morlet Zscore Power vs ' cSTATE_LBL{iBnd}]);
        
        % Bipolar Derivation
    subplot(2, 3, 6);
    CBASS_Plot_LFP_MorletPower(diff(db2LFP_Tst), sREC.inSampleRate/inFactor, false, bl1State_DS);
    title([cBAND_LABEL{iBnd} 'Morlet BiPol Power vs ' cSTATE_LBL{iBnd}]);
    
    CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
    close all;
    
end
% end