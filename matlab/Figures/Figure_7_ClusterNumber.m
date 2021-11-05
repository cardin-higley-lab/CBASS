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
chOutPath   = fullfile(chOutDir, 'Figure_7_ClusterNumber');

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDataPath);
% sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Sets general parameter
chLabel         = 'Gamma';
db1Band         = [30 80];
bl1State        = sREC.bl1Run;
chStateLabel    = 'Running';

% Sets options for enriched region identification
blUseRate   = false;
inMethod    = 1;

% Sets the number of cluster to be tested
in1NClu     = [5 20 100 200];
inNCnd      = length(in1NClu);

% Sets options for trough extraction
chDataFormat    = 'complex';
blZScore        = true;
inRefChan       = 5;

% Sets the gloabal label for the experiment;
chExpLbl    = [cEXP{iExp} '_' chLabel];
%% Peferms computation
% Gets the troughs
sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Initializes aggregation variables
db1Threshold        = nan(size(in1NClu));
[cSCORE, cEVENT]    = deal(cell(size(in1NClu)));

% Initializes the first figure
clear hFIG cFIG_NAME
hFIG(1) = figure('Position', [50 50 900 500]);
cFIG_NAME{1} = [chExpLbl '_CSD_Score_Overlap'];

% Loops through the number of cluster
for iCnd = 1:inNCnd
    % Sets the number of cluster for k-means partitioning
    inNClu  = in1NClu(iCnd);
    
    %Finds partitions
    fprintf('%d Events: ');
    tic
    sPART  = CBASS_L2_PartitionTrough(sTROUGH, bl1State, ...
        ~bl1State, blZScore, inNClu);
    toc
    
    % Computes a boolean indexing event (the final selection of troughs)
    bl1Event    = false(1, size(sREC.db2LFP, 2));
    bl1Event(sTROUGH.in1Index(sPART.bl1Partition)) = true;
    
    % Aggregate filter number and pulse traces
    cSCORE{iCnd}        = sort(sPART.db1Score);
    db1Threshold(iCnd)  = sPART.dbThreshold;
    cEVENT{iCnd}        = find(bl1Event);
    
    % Plots the CSD around pulse
    subplot(2, 4,  iCnd);
    CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, cEVENT{iCnd}, [-1 1] ./ mean(db1Band), sREC.inSampleRate);
    if iCnd == 1, title(sprintf('CSD around Events: %d Cluster', inNClu));
    else, title(sprintf('%d Clusters', inNClu)); end
end

%Calculates the event overlap 
[db2D_Asym, db2D_Sym] = EventOverlap(cEVENT, .005, sREC.inSampleRate);
%% Initializes the second figure

% Plost scores as a function of cluster number
subplot(2, 4, 5:6), hold on
for iCnd = 1:inNCnd
    db1Color    = [.2 1 .5] * iCnd / inNCnd;
    hPLT(iCnd)  = plot(cSCORE{iCnd}, 'DisplayName', num2str(in1NClu(iCnd)), ...
        'Color', db1Color);
    dbY = db1Threshold(iCnd);
    dbX = find(cSCORE{iCnd} > dbY, 1, 'first');
    plot([1 dbX], [dbY dbY], '--', 'Color', db1Color);
    plot([dbX dbX], [0 dbY], '--', 'Color', db1Color);
end
xlabel('Candiate Events'); ylabel('Probability Score');
legend(hPLT, 'Location', 'NorthWest');
title('Probability Score and Threshold vs Cluster Number');

% Plots the overlap between pulses
subplot(2, 4, 7:8)
imagesc(db2D_Asym);
set(gca, 'XTick', 1:inNCnd, 'XTickLabel', in1NClu, ...
    'YTick', 1:inNCnd, 'YTickLabel', in1NClu);
xlabel('N cluster'); ylabel('N cluster');
colorbar; title('Event Overlap');
%% Saves the figures
CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME, 'png');
close all;