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
chOutPath   = fullfile(chOutDir, 'Figure_7_ClusterNumber');

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

% Sets the number of cluster to be tested
in1NClu     = [25 50 100 200];
inNCnd      = length(in1NClu);

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

% Initializes aggregation variables
in1NFlt = nan(size(in1NClu));
cPULSE  = cell(size(in1NClu));

% Initializes the first figure
clear hFIG cFIG_NAME
hFIG(1) = figure('Position', [50 50 1250 600]);
cFIG_NAME{1} = [chExpLbl '_Filter1_CSD'];

% Loops through the number of cluster
for iCnd = 1:inNCnd
    % Sets the number of cluster for k-means partitioning
    inNClu  = in1NClu(iCnd);
    
    % Identifies enriched regions
    sFILTER = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, bl1State, blZScore, 1, inNClu);
    
    % Get the pulses for each enriched region
    sPULSE  = CBASS_L3_GetPulse_2(sREC, sFILTER(1).db2Filter, true, true);
    
    % Aggregate filter number and pulse traces
    in1NFlt(iCnd)   = length(sFILTER);
    cPULSE{iCnd}    = find(sPULSE(1).bl1Pulse);
    
    % Plots the filter
    subplot(2, inNCnd, iCnd)
    imagesc(sFILTER(1).db2Filter);
    title(sprintf('Filter: %d Cluster', inNClu));
    
    % Plots the CSD around pulse
    subplot(2, inNCnd, inNCnd + iCnd);
    CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sPULSE(1).bl1Pulse, [-1 1] ./ mean(db1Band), sREC.inSampleRate);
    title(sprintf('CSD around pulse: %d Cluster', inNClu));
end

% Initializes the first figure
hFIG(2) = figure('Position', [50 50 1250 600]);
cFIG_NAME{2} = [chExpLbl '_NFilt_Overlap'];

% Plots the number of filter as a function of the number of k-means cluster
subplot(1, 2, 1);
plot(1:inNCnd, in1NFlt, '--ok')
set(gca, 'XTick', 1:inNCnd, 'XTickLabel', in1NClu);
xlim([0 inNCnd + 1]);
xlabel('K-Means cluster'); ylabel('Number of filters');

% Plots the overlap between pulses
[db2D_Asym, db2D_Sym] = EventOverlap(cPULSE, .005, sREC.inSampleRate);
subplot(1, 2, 2);
imagesc(db2D_Asym);
set(gca, 'XTick', 1:inNCnd, 'XTickLabel', in1NClu, ...
    'YTick', 1:inNCnd, 'YTickLabel', in1NClu);
xlabel('K-Means cluster'); ylabel('K-Means cluster');
colorbar; title('Event Overlap');

% Saves the figures
CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
close all;