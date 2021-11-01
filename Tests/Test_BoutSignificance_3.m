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
chOutPath   = fullfile(chPath,'Tests', 'Test_BoutSignificance_3');
chTmpPath   = fullfile(chPath,'Tests', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
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

% Sets options for the phase randomized signal and the normalization of
% pulse detection
in1Norm     = [1 2];
cCNTR_TYPE  = {'Min', 'Mean'};
inNCnd      = length(in1Norm);
%% Loops through conditons
for iBnd = 1:length(cBAND)
    % Sets band labels
    chLabel     = cBAND_LABEL{iBnd};
    db1Band     = cBAND{iBnd};
    
    % Sets the gloabal label for the experiment;
    chExpLbl    = [chExperiment '_' chLabel];
    
    % Initializizes the figure
    clear hFIG
    hFIG(1) = figure('Position', [50 50 1250 600]);
    hFIG(2) = figure('Position', [50 50 1250 600]);
    hFIG(3) = figure('Position', [50 50 1250 600]);
%     hFIG(4) = figure('Position', [50 50 1250 600]);
    cFIG_NAME = {[chExpLbl '_PeakHistogram'], ...
        [chExpLbl '_CSDAroundPulse'], [chExpLbl '_PulseEnrichment']};
    
    % Gets the troughs
    sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
    
    % Gets the rand trough
    sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
    
    % Identifies enriched regions
    sFILTER     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, cSTATE{iBnd}, blZScore, 1, inNClu);
    if isempty(sFILTER), sprintf('No Filter for %s %s Method %s', chExperiment, cBAND_LABEL{iBnd}, iMth); continue; end
    
    % Loops through condition to define pulse
    for iCnd = in1Norm
        % Sets the loop condition label.
        chCondLbl   = sprintf('%s_Centering', cCNTR_TYPE{iCnd});
        
        % Get the pulses for each enriched region
        if iCnd == 1, sPULSE = CBASS_L3_GetPulse_3(sREC, sFILTER(1).db2Filter, true, true);
        elseif iCnd == 2, sPULSE = CBASS_L3_GetPulse_2(sREC, sFILTER(1).db2Filter, true, false); end
        
        % Plots the pulse value histogram
        figure(hFIG(1)), subplot(1, 2, iCnd);
        CBASS_Plot_PeakHistogram(sPULSE);
        ylabel(strrep(chCondLbl, '_', ' '))
        
        % Plots the pulse triggerred averages
        figure(hFIG(2)), subplot(1, 2, iCnd);
        CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sPULSE.bl1Pulse, [-2 2]./mean(db1Band), sREC.inSampleRate);
        title(strrep(chCondLbl, '_', ' '))
        
        % Plots event enrichment
        figure(hFIG(3)), subplot(1, 2, iCnd)
        chSigString = CBASS_Plot_Pulse_vs_State(sPULSE.bl1Pulse, cSTATE{iBnd}, [-5 10], sREC.inSampleRate);
        title(sprintf('%s : %s', strrep(chCondLbl, '_', ' '), chSigString))
    end
    CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
    close all;
end
end