%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
for iExp = 1:length(cEXP)
% iExp = 1
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Figures_and_Tests', 'Test_Pipeline');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDir);
sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Sets general parameter
% Sets options for trough extraction
chDataFormat    = 'complex';
blZScore        = true;
inRefChan       = 5;

% Sets options for enriched region identification
blUseRate   = false;
inMethod    = 1;
inNClu      = 20;

%% Sets the loop
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND       = {[15 30], [30 80]};
cSTATE      = {sREC.bl1Pres sREC.bl1Run};
cSTATE_LBL  = {'Stim', 'Running'};

% Loops through conditons
for iBnd = 1:length(cBAND)
    % Sets band labels
    chLabel     = cBAND_LABEL{iBnd};
    db1Band     = cBAND{iBnd};
    
    % Sets the gloabal label for the experiment;
    chExpLbl    = [chExperiment '_' chLabel];
    
    % Gets the troughs
    sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
    
    % Gets the rand trough
    sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
    
    % Identifies enriched regions
    sFILTER     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, cSTATE{iBnd}, blZScore, 1, inNClu);
    if isempty(sFILTER), sprintf('No Filter for %s %s Method %s', chExperiment, cBAND_LABEL{iBnd}, 1); continue; end
    
    % Get the pulses for each enriched region
    sPULSE = CBASS_L3_GetPulse_2(sREC, sFILTER(1).db2Filter, true, true);

    % Initializizes the figure
    clear hFIG
    hFIG = figure('Position', [50 50 1250 600]);
    cFIG_NAME = {chExpLbl};
    
    % Plots the spectra
    subplot(2, 2, 1);
    CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate, false, cSTATE{iBnd});
    title([cBAND_LABEL{iBnd} ' Power vs ' cSTATE_LBL{iBnd}]);
    
    % Plots the pulse triggerred averages
    subplot(2, 2, 2);
    CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sPULSE.bl1Pulse, [-2 2]./mean(db1Band), sREC.inSampleRate);
    title([strrep(chExpLbl, '_', '\_') ' CSD around pulse']);
    
    % Plots event enrichment
     subplot(2, 2, 3)
    chSigString = CBASS_Plot_Pulse_vs_State(sPULSE.bl1Pulse, cSTATE{iBnd}, [-5 10], sREC.inSampleRate);
    title(sprintf('%s', chSigString))
    
    % Plots the filters as a function under and above the 80%percentile of
    % a moving average of the number of pulses
    dbWinLenSec     = 0.5;
    db1Win          = rectwin(dbWinLenSec * sREC.inSampleRate);
    db1NPulse       = conv(sPULSE.bl1Pulse, db1Win, 'same')./dbWinLenSec;
    db80pct         = quantile(db1NPulse, .8);
    bl180pctPlus    = db1NPulse > db80pct;
    subplot(2, 2, 4);
    CBASS_Plot_LFP_FourierPower(sREC.db2LFP, sREC.inSampleRate,...
        false, bl180pctPlus, {sprintf('<%.1fHz', db80pct) sprintf('>%.1fHz', db80pct)});
    title([cBAND_LABEL{iBnd} ' Power vs Pulse Frequency'])
    
    CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
    close all;
end
end