%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
for iExp = 1:length(cEXP)
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Tests', 'Test_Beta_Gamma_Overlap');
chTmpPath   = fullfile(chPath,'Tests', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

%% Sets data formatting and embeding global parameters
% Sets the data format
chDataFormat    = 'complex';

% Call python to do the embeding
blZScore        = true;
if blZScore, chZFlag = 'ZScore'; else, chZFlag = 'Raw'; end

% Sets the embeding of choice for
chEmbedMethod   = 'umap';
inN_Component   = 3;
%% Loads the data
sREC        = CBASS_L0_LoadData(chDir);
sREC        = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Set up the loop
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND       = {[15 30], [30 80]};
cSTATE      = {sREC.bl1Pres sREC.bl1Run};
cSTATE_LBL  = {'Stim', 'Running'};
inNBnd      = length(cBAND_LABEL);

% Sets option
blUseRate   = false;
inMethod    = 1;
inNClu      = 20;
%% Loops through conditons
% Initizalize aggregate variable
cEVENT = {};

for iBnd = 1:inNBnd
    % for iBnd = 2
    chLabel         = cBAND_LABEL{iBnd};
    dbBand          = cBAND{iBnd};
    inRefChan       = 5;
    
    % Sets the gloabal label for the experiment
    chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat '_' chZFlag];
    chExpLblShrt    = [chExperiment '_' chLabel];
    
    % START PROCESSING ------------------------------------------------
    % Computes the trough
    sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
    
    % Identifies enriched regions
    [sREGION, sCLU] = CBASS_L2_GetEnrichedRegion(sTROUGH, cSTATE{iBnd}, blZScore, blUseRate, inMethod, inNClu);
    if isempty(sREGION); fprintf('No enriched region for %s\r', cBAND_LABEL{iBnd}); continue; end
    
    % Get the pulses for each enriched region
    sPULSE = CBASS_L3_GetPulse(sREC, sTROUGH, sREGION(1).bl1Member);
    
    % Aggregates the pulse events with the most enriched region
    cEVENT = cat(2, cEVENT, {sPULSE.bl1Pulse});
end

% 1 --- Plots overlap image
hFIG = figure('Position', [50 50 1250 600]);
chFigName = ['Overlap_Beta_Gamma_' chExpLblShrt];
cEVENT_IDX = cellfun(@find, cEVENT, 'UniformOutput', false);
db2D_Asym = EventOverlap(cEVENT_IDX, .005, sREC.inSampleRate);
imagesc(1:inNBnd, 1:inNBnd, db2D_Asym)
ylabel('Cond'); xlabel('Cond'); title('EventOverlap');
set(gca, 'XTick', 1:inNBnd, 'XTickLabel', cBAND_LABEL, 'YTick', 1:inNBnd, 'YTickLabel', cBAND_LABEL);
xtickangle(90);
colorbar
CBASS_SaveFig(chOutPath, hFIG, {chFigName});

% 2 -- Plots LFP example around stim onset
in1Anchor       = find(cSTATE{1}(1:end - 1) & cSTATE{1}(2:end));

hFIG = figure('Position', [50 50 1250 600]);
chFigName = ['LFP_Examples_BetaGamma_' chExpLabel];
for iEx = 1:4
    subplot(2, 2, iEx)
    CBASS_Plot_LFP_EventExample(sREC.db2LFP, in1Anchor(iEx), [-2 3], sREC.inSampleRate, cEVENT, cBAND_LABEL)
    title(sprintf('Data around %s onset: %d', cSTATE_LBL{1}, iEx));
end
CBASS_SaveFig(chOutPath, hFIG, {chFigName});

close all
end