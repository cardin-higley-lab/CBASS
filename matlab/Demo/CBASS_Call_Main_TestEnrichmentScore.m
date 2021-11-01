clear, close all
%% Adds the code to the path ---------- change to the path to the code on your computer -----------------
chCodePath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/Pipeline'; %'D:\gamma_bouts\Pipeline';
chUtilPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/Pipeline/Utilities'; %'D:\Utilities';
addpath(genpath(chCodePath)); %Add the root folder to current folder
addpath(genpath(chUtilPath)); %Add utilities to the path
%% Load the data ----- change to load your data ---------------------------------------------------------
% The section must output the variable db2LFP (i.e.  a channel x
% time_sample matrix) and its sample rate inSample rate

% Path to the data
chDataPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/Data/BCP_CBASS_Test';
chExperiment    = '2017-03-30_16-02-10_2';
chExpPath       = fullfile(chDataPath, chExperiment);
chMDSFile       = 'BCP_L4_MakeMetaDataStructure.mat';

%Load the input
sINPUT = load(fullfile(chExpPath, chMDSFile), '-mat'); sINPUT = sINPUT.sCFG;

%Get the overal sample rate
inSampleRate = sINPUT.sL4MMDS.inWorkSampleRate;

%Extract the LFP
db2LFP = sINPUT.sL4MMDS.db2LFP; %Exctracts the LFP matrix
db2LFP = PCP_U_NormalizeLFP(db2LFP, inSampleRate); % Filters and z-scores the LFP
db2LFP = db2LFP(2:end, :); % Removes the first channel of the LFP which is set as the reference. We are left with 15 channels.

%Get the indices of when the mouse is running or whisking
bl1Run          = sINPUT.sL4MMDS.bl1WheelOn;

%Extracts the stimulation times
bl1FullCtr      = [sINPUT.sL4MMDS.sTRIAL.dbContrast] == 100;
in1PresOnIdx    = NS_GetTStampEventIndex(sINPUT.sL4MMDS.db1TStamps, [sINPUT.sL4MMDS.sTRIAL(bl1FullCtr).dbStimOnsetTS]);
in1PresOffIdx   = NS_GetTStampEventIndex(sINPUT.sL4MMDS.db1TStamps, [sINPUT.sL4MMDS.sTRIAL(bl1FullCtr).dbStimOffsetTS]);
bl1Stim         = NS_MakeEpochVector(in1PresOnIdx, in1PresOffIdx, length(sINPUT.sL4MMDS.db1TStamps));

%Extracts lick times
bl1LickBout 	= sINPUT.sL4MMDS.bl1LickBout;
bl1LickExcl		= sINPUT.sL4MMDS.bl1LickExcl;

%Computes states
bl1RunOnly 		= bl1Run & ~bl1Stim & ~bl1LickBout & ~bl1LickExcl;
% bl1StimOnly 	= ~bl1Run & bl1Stim & ~bl1LickBout & ~bl1LickExcl;
% bl1LickOnly   = ~bl1Run & ~bl1Stim & bl1LickBout & ~bl1LickExcl;
bl1Quiet		= ~bl1Run & ~bl1Stim & ~bl1LickBout & ~bl1LickExcl;

sREC.db2LFP         = db2LFP;
sREC.inSampleRate   = inSampleRate;
%% Sets the output folder ----- edit or comment to provide externally throug a system call -------------
chOutFile   = 'Test_Enrichment score';
chOutPath   = fullfile('D:\BCP_CBASS_Test\', chOutFile);
if ~exist(chOutPath, 'dir'), mkdir(chOutPath); end
%% Non optional parameters ----- edit or comment to provide externally throug a system call -------------
cBAND           = [30 80];
bl1State        = bl1RunOnly;
%% Optional parameters ----- Uncomment to provide externally through a system call using carrier name --
sOPTION.cBAND_LABEL     = {'Gamma'};  % character array or cell array of labels  for bands in cBAND (i.e. {'Beta', 'Gamma'})
sOPTION.cSTATE_LBL      = {'Running'};  % character array or cell array of labels  for states in cSTATE (i.e. {'Stim', 'Running'})
sOPTION.blVerbose       = true;    % Sets whether to print progress on the screen. Default is true

%L1 options for formatting hilbert troughs (Function CBASS_L1_GetTrough)
sOPTION.chDataFormat    = 'complex'; % Format of hilbert transform coordinates ('polar' or 'complex' Default is 'complex')
sOPTION.inRefChan       = 5;    % Reference channel (Default is the last of row of db2LFP.

%L2 options for the computation of filters (Function CBASS_L2_GetFilters)
sOPTION.cBASELINE       = bl1Quiet; % Baseline state.
sOPTION.blZScore        = true;     % z-score for k-means partitioning. Default is true.
sOPTION.inMethod        = 1;     % method for spectral clusering of template. Default is 1.
sOPTION.dbSigThrs       = 10.^-4;    % threshold for enrichment significance. Default is 10.^-4.
sOPTION.inNMaxIter      = 10000;   % maximum iteration of k-means partioning. Default is 10000;

%L3 options for template matching
sOPTION.blCntr          = true;       % Sets whether to center signal and template about their mean. Default is true
sOPTION.blNorm          = true;       % Sets whether to normalize signal and template by their S.D. Default is true.
%% Calls main function
% if ~exist('sOPTION', 'var'), sOPTION = struct; end % Creates sOPTION if
% it does not exist sFREQ_BAND =
% CBASS_Main_DetectEnrichedMotif(sREC.db2LFP, sREC.inSampleRate, cBAND,
% cSTATE, sOPTION);
%% Calls the beginning of the pipeline
% Adds phase randomized signal
if sOPTION.blVerbose, fprintf('---->> Add phase randomized signal ... '); end
sREC                = CBASS_L1_AddPhaseRandomizedSignal(sREC);

% Extracts trougths for the real and surrogates signals
if sOPTION.blVerbose, fprintf('---->> Extract hilbert troughs ... '); end
sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, inSampleRate, cBAND, ...
    sOPTION.inRefChan, sOPTION.cBAND_LABEL, sOPTION.chDataFormat);
sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, inSampleRate, cBAND, ...
    sOPTION.inRefChan, sOPTION.cBAND_LABEL, sOPTION.chDataFormat);
if sOPTION.blVerbose, fprintf('Done\r\r'); end

% Initializes figure numbers
iFig = 0;
clear hFIG cFIG_NAME
%% Test the enrichment score algorithm
db2Trough   = sTROUGH.db2Trough;
bl1T_State  = bl1State(sTROUGH.in1Index);
dbSigThrs   = .05;
in1NItr     = [100 300 1000 3000];
[inRow, inCol] = FindPlotNum(length(in1NItr));

% 
% for iThr = [10, 20, 40, 80]
%     
%     iFig = iFig + 1;
%     hFIG(iFig) = figure('Position', [50 50 1250 600]);
%     cFIG_NAME{iFig} = sprintf('%d Cluster', iThr);
%     
%     for iPlt = 1:length(in1NItr)
%         inNItr = in1NItr(iPlt);
%         tic, in1Score = CBASS_U_EnrichmentScore(db2Trough, bl1T_State, ~bl1T_State, iThr, dbSigThrs, true, inNItr); toc
%         tic, in1Score2 = CBASS_U_EnrichmentScore(db2Trough, bl1T_State, ~bl1T_State, iThr, dbSigThrs, true, inNItr); toc
%         
%         [in1Srt_1, in1Idx] = sort(in1Score);
%         in1Srt_2 = in1Score2(in1Idx);
%         
%         subplot(inRow, inCol, iPlt)
%         plot(in1Srt_2); hold on, plot(in1Srt_1)
%         title(sprintf('%d Clu: %d iteration', iThr, inNItr));
%     end
% end
%% Test significance levels
% Sets data and data rand
db2Trough   = sTROUGH.db2Trough;
bl1T_State  = bl1State(sTROUGH.in1Index);
db2TrghRnd  = sTRGH_RND.db2Trough;
bl1T_StRnd  = bl1State(sTRGH_RND.in1Index);
dbSigThrs   = .05;

% Set number of interation and cluster number
inNClu = 20;
inNItr = 1000;

% Finds troughs situated in enriched regions (Preselct the trougths)
in1Score    = CBASS_U_EnrichmentScore(db2Trough, bl1T_State, ~bl1T_State, inNClu, dbSigThrs, true, inNItr);
in1ScrRnd   = CBASS_U_EnrichmentScore(db2TrghRnd, bl1T_StRnd, ~bl1T_StRnd, inNClu, dbSigThrs, true, inNItr);
dbSig_1     = quantile(in1ScrRnd, .95);
dbSig_2     = .95;
%%
% Computes a histogram of the rate of enrichment as a fonction of the
% sorted score
[in1SrtScore, in1SortIdx] = sort(in1Score);
bl1SrtT_State   = bl1T_State(in1SortIdx);
db1StateRate    = conv(bl1SrtT_State, rectwin(1000), 'same')/1000;
[in1SrtScrRnd, in1SortIdx] = sort(in1ScrRnd);
bl1SrtT_State   = bl1T_StRnd(in1SortIdx);
db1StRtRnd      = conv(bl1SrtT_State, rectwin(1000), 'same')/1000;

% % Initializes the figure
% iFig = iFig + 1;
% hFIG(iFig) = figure('Position', [50 50 1250 600]);
% cFIG_NAME{iFig} = 'SigLevels';
% 
% % Plots the enrichment score as well as the different levels of
% % significance
% subplot(1, 2, 1)
% plot(in1SrtScore); hold on
% plot(in1SrtScrRnd);
% plot(db1StRtRnd);
% legend({'SortedScore', 'SortedRndScore', 'RandRate'}, 'Location', 'NorthWest');
% title(sprintf('Enrichment Score: %d Clu, %d iteration', inNClu, inNItr))
% 
% % Plots the enrichment score as well as the different levels of
% % significance
% subplot(1, 2, 2)
% plot(in1SrtScore); hold on
% plot(db1StateRate);
% plot(xlim, [1 1] * dbSig_1, 'r--');
% plot(xlim, [1 1] * dbSig_2, 'g--');
% legend({'SortedScore', 'StateRate', '95pct on Rand', '95pct'}, 'Location', 'NorthWest');

%% Test the enrichment as a function fo the threshold
% Exctracts data
db2Trough       = sTROUGH.db2Trough;
bl1T_State      = bl1State(sTROUGH.in1Index);

% Sets parameters
dbSigThrs       = 0.05;
inNClu          = 20;
inNItr          = 1000;

% Computes enrichment score
in1Score        = CBASS_U_EnrichmentScore(db2Trough, bl1T_State, ~bl1T_State, inNClu, dbSigThrs, true, inNItr);

% Sets the number of clusters
db1Thrs = [.95 .90 .80];
inNThrs = length(db1Thrs);

% Sets the number of plots
[inNRow, inNCol] = FindPlotNum(inNThrs + 1);

% Initializes the figure
iFig            = iFig + 1;
hFIG(iFig)      = figure('Position', [50 50 1250 600]);
cFIG_NAME{iFig} = 'Rate_vs_Threshold';

% Loops through cluster numbers
for iThr = 1:length(db1Thrs)
    
    % Plots the enrichement of the trought identified in the filter
    subplot(inNRow, inNCol, iThr)
    bl1Sel = in1Score > db1Thrs(iThr);
    chSigString = CBASS_Plot_Pulse_vs_State(sTROUGH.in1Index(bl1Sel), bl1Run, [-5 10], sREC.inSampleRate);
    title(sprintf('Thrs: %.2f: %s', db1Thrs(iThr), chSigString))
end

% Plots the natural rate of events
% Plots the enrichement of the trought identified in the filter
subplot(inNRow, inNCol, iThr + 1)
dbThrsNat   = quantile(in1Score, 1 - (sum(in1Score)/length(in1Score)));
bl1Sel      = in1Score > dbThrsNat;
chSigString = CBASS_Plot_Pulse_vs_State(sTROUGH.in1Index(bl1Sel), bl1Run, [-5 10], sREC.inSampleRate);
title(sprintf('Natural Rate - Thrs: %.2f: %s', dbThrsNat, chSigString));

%% Plots the proportion of events retained as a function of the threhosld

% Initializes the figure
iFig            = iFig + 1;
hFIG(iFig)      = figure('Position', [50 50 1250 600]);
cFIG_NAME{iFig} = 'Event_Metrics_vs_Threshold';

% Plost the scores
subplot(2, 2, 1); plot(sort(in1Score));
title('Sorted enrichment score'); 
ylabel('Score (fraction of trials where enriched');
xlabel('Event (i.e. trougths)');

% Plots the fraction of events retained as a function of threshold
subplot(2, 2, 2); plot(sort(in1Score), (length(in1Score):-1:1)/length(in1Score));
hold on;
plot([1 1] * dbThrsNat, [0 1],  '--r');
for iThr = 1:inNThrs, plot([1 1] * db1Thrs(iThr), [0 1], '--g'); end
xlabel('Threshold'), ylabel('Fraction of events retained')

% Plots the fraction of events retained as a function of threshold
subplot(2, 2, 3);
db1Thr_X        = 0:.01:.99;
db1ONonOFF_Y    = nan(size(db1Thr_X)); 
for iThr = 1:length(db1Thr_X)
    bl1Sel = in1Score > db1Thr_X(iThr);
    db1ONonOFF_Y(iThr) = (sum(bl1T_State(bl1Sel)) * inSampleRate /sum(bl1State))  ...
        /(sum(~bl1T_State(bl1Sel)) * inSampleRate /sum(~bl1State));
end
plot(db1Thr_X, db1ONonOFF_Y); hold on
plot([1 1] * dbThrsNat, ylim,  '--r');
for iThr = 1:inNThrs, plot([1 1] * db1Thrs(iThr), ylim, '--g'); end
xlabel('Threshold'), ylabel('Enrichment (Rate during vs outside run)');

% Plots the fraction of events retained as a function of threshold
subplot(2, 2, 4); plot(sort(in1Score), (length(in1Score):-1:1) * inSampleRate /size(db2LFP, 2));
hold on;
plot([1 1] * dbThrsNat, ylim,  '--r');
for iThr = 1:inNThrs, plot([1 1] * db1Thrs(iThr), ylim, '--g'); end
xlabel('Threshold'), ylabel('Event Rate (Hz)')
%% Sav/es the figures
CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);