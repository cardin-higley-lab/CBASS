%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};

iExp = 1;
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'Data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath, 'Figures', 'Figure_5_6_TemplateClustering');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDir);
sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);

% Sets general parameter
chLabel         = 'Beta';
db1Band         = [15 30];
bl1State        = sREC.bl1Pres;
chStateLabel    = 'Stim';

% Sets options for enriched region identification
blUseRate   = false;
inMethod    = 1;
inNClu      = 20;

% Sets options for trough extraction
chDataFormat    = 'complex';
blZScore        = true;
inRefChan       = 5;

% Sets the gloabal label for the experiment;
chExpLbl    = [chExperiment '_' chLabel];

% Gets the troughs
sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Gets the rand trough
sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Identifies enriched regions
[sFILTER, sCLU]     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, bl1State, blZScore, 2, inNClu);

% Makes filters for significantly enriched regions
%Finds signigicantly enriched regions
in1Sel              = find([sCLU.dbRate_Dev] > 0 & [sCLU.dbPVal] < 10.^-4);
dbCycLen            = 1./mean(db1Band);
db1Filter_WinSec    = [-dbCycLen dbCycLen] / 2;

% Get the filter for each enriched region
cFILTER = {}; db2FiltMat = [];
for iClu = in1Sel
    % Makes the filter
    in1EventIdx = sTROUGH.in1Index(sCLU(iClu).bl1Member);
    db2Filter   = CBASS_U_MakeFilter(sREC.db2LFP, sREC.inSampleRate, in1EventIdx, db1Filter_WinSec);
    
    % Aggregates the pulse events for the cluster of interest
    cFILTER     = cat(2, cFILTER, {db2Filter});
    db2FiltMat  = cat(2, db2FiltMat, db2Filter(:));
end

% Plots the first figure
% Initializizes the figure
clear hFIG cFIG_NAME
hFIG(1) = figure('Position', [50 50 1250 600]);
cFIG_NAME{1} = [chExpLbl '_5_EventTriggeredAverages'];

% Plots the figure
inNFlt  = length(cFILTER);
[inNRow, inNCol] = FindPlotNum(inNFlt);
for iFlt = 1:inNFlt
    subplot(inNRow, inNCol, iFlt)
    imagesc(cFILTER{iFlt}); title(sprintf('Filter %d', iFlt));
end

%Calculate the correlation between fitlters and sets those inferior
%to the threshold to zero to built the adjacency matrix
db2Adj  = corr(db2FiltMat);

% Uses MDS to represent filters according to their coordinates
db1XY   = cmdscale(db2Adj, 2);

% Initializizes the figure
hFIG(2) = figure('Position', [50 50 1250 600]);
cFIG_NAME{2} = [chExpLbl '_SpectralClustering'];

% Plots the fuly connected graph
subplot(1, 3, 1); hold on
plot(db1XY(:, 1), db1XY(:, 2), 'kx')
for iFlt1 = 1:inNFlt - 1
    for iFlt2 = iFlt1 + 1:inNFlt
        plot(db1XY([iFlt1 iFlt2], 1), db1XY([iFlt1 iFlt2], 2), 'Color', [.8 .8 .8]);
    end
end
for iFlt = 1:inNFlt
    text(db1XY(iFlt, 1), db1XY(iFlt, 2), ['  ' num2str(iFlt)])
end
title('Fully Connected Graph');

% Plots the number of cluster as a function of the threshold between the
% graphs
db1Thrs = .01:.01:.99;
in1NFlt = zeros(size(db1Thrs));
for iThr = 1:length(db1Thrs)
    % Trims edges as a function of their threshold
    db2Adj_Thr = db2Adj;
    db2Adj_Thr(db2Adj_Thr < db1Thrs(iThr)) = 0;
     
    % Gets the number of clusters for spectral clustering
    in1CluSC        = CBASS_U_SpectralCluster(db2Adj_Thr);
    in1NFlt(iThr)   = max(in1CluSC);
end
%Calculates a threshold for the purpose of illustration (the real threshold
%would yield one cluster)
inNFlt_Fig      = 3;
dbThreshold     = mean(db1Thrs(in1NFlt == inNFlt_Fig));  % Threshold
db2Adj_Thr = db2Adj;
db2Adj_Thr(db2Adj_Thr < dbThreshold) = 0;
in1Clu          = CBASS_U_SpectralCluster(db2Adj_Thr);   % Cluser attribution at that treshold
    % Does the actual plot
subplot(1, 3, 2); hold on
db1YL = [0 in1NFlt(end) + 1];
plot(db1Thrs, in1NFlt); hold on;
plot([1 1] * dbThreshold, db1YL, '--r');
xlim(db1Thrs([1 end])); ylim(db1YL);
xlabel('Similarity threshold'); ylabel('Number of filters');

% Plots the subgraphs
subplot(1, 3, 3); hold on
plot(db1XY(:, 1), db1XY(:, 2), 'kx')
for iFlt = 1:inNFlt_Fig
    in1CluSel   = find(ismember(in1Clu, iFlt));
    inNFltSel   = length(in1CluSel);
    for iFlt1 = 1:inNFltSel - 1
        for iFlt2 = iFlt1 + 1:inNFltSel
            plot(db1XY(in1CluSel([iFlt1 iFlt2]), 1), db1XY(in1CluSel([iFlt1 iFlt2]), 2), 'Color', [.8 .8 .8]);
        end
    end
end
for iFlt = 1:inNFlt
    text(db1XY(iFlt, 1), db1XY(iFlt, 2), ['  ' num2str(iFlt)])
end
title('Trimmed Graph');

CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
close all;