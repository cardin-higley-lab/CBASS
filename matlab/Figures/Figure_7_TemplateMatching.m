%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};

iExp = 2;
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath, 'Figures', 'Figure_7_TemplateMatching');
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
sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Gets the rand trough
sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

% Identifies enriched regions
sFILTER     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, bl1State, blZScore, 1, inNClu);

% Get the pulses for each enriched region
sPULSE      = CBASS_L3_GetPulse_2(sREC, sFILTER(1).db2Filter, true, true);

% Initializizes the figure
clear hFIG
hFIG = figure('Position', [50 50 1250 600]);
cFIG_NAME = {chExpLbl};

% Plots the pulse triggerred averages
in1StimON   = find(~bl1State(1:end - 1) & bl1State(2:end));
inAnchor    = in1StimON(2) - (0.2 * sREC.inSampleRate);
db1WinSec   = [0 .6];
subplot(1, 2, 1);
CBASS_Plot_LFP_EventExample(sREC.db2LFP, inAnchor, db1WinSec, sREC.inSampleRate);
% Appends the filter
db1YL       = ylim;
db1XL       = xlim;
inNChan     = size(sREC.db2LFP, 1);
dbScale     = db1YL(2) ./ (inNChan + 1.5);
db2Flt      = sFILTER(1).db2Filter;
inNSmpFlt   = size(db2Flt, 2);
db1FTime    = (-2 * inNSmpFlt: -1 * inNSmpFlt - 1) ./ sREC.inSampleRate;
% Plots the LFP
hold on
for iChan = 1:inNChan
    plot(db1FTime, (dbScale * (inNChan - iChan + 1)) + 2 * db2Flt(iChan, :), 'r');
end
title('Filter(Red) and LFP (black)');
xlim([db1FTime(1) db1XL(2)])


% Sets parameters for plotting score and pulse exerpt
in1SelIdx   = inAnchor + (round(db1WinSec(1) * sREC.inSampleRate):round(db1WinSec(2) * sREC.inSampleRate));
db1X        = (in1SelIdx - in1SelIdx(1) - 1) ./ sREC.inSampleRate;
db1ScoreSel = sPULSE(1).db1Score(in1SelIdx);
in1PeakIdx  = find(sPULSE(1).bl1Pulse(in1SelIdx));
% Plots the score it's envelope and the peaks
subplot(1, 2, 2);
plot(db1X, db1ScoreSel); hold on
plot(db1X(in1PeakIdx), db1ScoreSel(in1PeakIdx), 'rx')
plot(db1X([1 end]), sPULSE.dbThrPeak * [1 1], '--')
title('Template Matching'); ylabel('Score Exerpt (A.U.)'); xlabel('Time (s)');
legend({'Score', 'Pulse', 'Threshold'});
ylim([-.7 1]); xlim(db1X([1 end]));

% Saves figures
CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
close all;