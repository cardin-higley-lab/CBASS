% Runs the data analysis script
cd 'D:\gamma_bouts\scripts'; % < -------------Edit
run LoadDataScript_BetaGamma_HilbertFormatting.m
close all
load('ClusterBranchIndices_v2.mat');

% Creates the figure directory
chFigDir = 'BoutSignificance';
if ~exist(chFigDir, 'dir'); mkdir(chFigDir); end
iFig = 0;
clear hFIG
%%
% Defines Bands
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND = {[15 30], [30 80]};

% Defines what band we are dealing with here
iBnd = 2;

% Defines whether to center and normalize
blCntr = true;
blNorm = true;

% Gets the size of the LFP
[inNChan, inNSample] = size(db2LFP);

% Gets the indices of the branches
in1BrchIdx  = in1BranchIdx_10Clusters;

% Gets the indices of the epoches on the initial LFP recording
in1TroughIdx    = cTROUGH_IDX{2};
in1EventIdx     = in1TroughIdx(in1BrchIdx == 1); % Clarifies event idx
inNEvt          = length(in1EventIdx);

% Defines the window
dbCycLen            = 1./mean(cBAND{iBnd});
db1Filter_WinSec    = [-dbCycLen dbCycLen] / 2;
in1Filt_RelIdx      = round(inSampleRate * db1Filter_WinSec(1)) : ...
    round(inSampleRate .* db1Filter_WinSec(2));

% Creates the ETA
db3LFP_ETA = nan(inNChan, length(in1Filt_RelIdx), inNEvt);
for iEvt = 1:inNEvt
    db3LFP_ETA(:, :, iEvt) = db2LFP(:, in1EventIdx(iEvt) + in1Filt_RelIdx);
end
db2Filter = nanmean(db3LFP_ETA, 3);

% Runs template matching, detects the peaks and computes the enveloppe
db1Score    = CBASS_U_MultiChannelTemplateMatching(db2LFP, db2Filter, blCntr, blNorm);
bl1Peak     = [false, diff(db1Score(1:end - 1)) > 0 & diff(db1Score(2:end)) < 0, false];
db1Env      = abs(hilbert(db1Score));

% Phase randomizes the LFP;
[db2Coeff, db2LFP_Rnd] = pca(db2LFP');
db2LFP_Rnd = CBASS_U_PhaseRandomize1D(db2LFP_Rnd, 1);
db2LFP_Rnd = (db2LFP_Rnd/db2Coeff)';

% Runs the template matching on the phase randomized LFP, detects the peaks
% and computes the enveloppe
db1Score_Rnd    = CBASS_U_MultiChannelTemplateMatching(db2LFP_Rnd, db2Filter, blCntr, blNorm);
bl1Peak_Rnd     = [false, diff(db1Score_Rnd(1:end - 1)) > 0 & diff(db1Score_Rnd(2:end)) < 0, false];
db1Env_Rnd      = abs(hilbert(db1Score_Rnd));
%% A Quick plot of the results
% Computes time selection for plot.
dbTimeSec   = 2;
dbTOffset   = 1215;
in1T_Idx    = (1:round(inSampleRate * dbTimeSec)) + round(inSampleRate * dbTOffset);
db1Time     = in1T_Idx ./ inSampleRate;

% Defines the bins for the histogram
if blNorm, db1Bin = -.5:.01:1; else; [~, db1Bin] = hist([db1Score(bl1Peak) db1Score_Rnd(bl1Peak_Rnd)], 150); end

% Creates the figure
iFig = iFig + 1; hFIG(iFig) = figure('Position', [50 50 1250 600]);
subplot(3, 2, 1)
db1Scr_Sel  = db1Score(in1T_Idx);
plot(db1Time, db1Scr_Sel), hold on
plot(db1Time, db1Env(in1T_Idx));
plot(db1Time(bl1Peak(in1T_Idx)), db1Scr_Sel(bl1Peak(in1T_Idx)), 'xr')
title('LFP'), xlabel('Time (s)'), ylabel('TM Score')
subplot(3, 2, 2)
db1Scr_Sel  = db1Score_Rnd(in1T_Idx);
plot(db1Time, db1Scr_Sel), hold on
plot(db1Time, db1Env_Rnd(in1T_Idx));
plot(db1Time(bl1Peak_Rnd(in1T_Idx)), db1Scr_Sel(bl1Peak_Rnd(in1T_Idx)), 'xr')
title('Phase Randomized LFP'), xlabel('Time (s)'), ylabel('TM Score')
subplot(3, 2, 3)
db1HistBout = hist(db1Score(bl1Peak), db1Bin);
bar(db1Bin, db1HistBout);
xlabel('LFP TM Score at Peak'); ylabel('Count')
subplot(3, 2, 4)
db1HistCBASS_Rnd = hist(db1Score_Rnd(bl1Peak_Rnd), db1Bin);
bar(db1Bin, db1HistCBASS_Rnd);
xlabel('Rand TM Score at Peak'); ylabel('Count')
subplot(3, 2, 5)
db1HistEnv = hist(db1Env, db1Bin);
bar(db1Bin, db1HistEnv);
xlabel('LFP TM Score Envelope'); ylabel('Count')
subplot(3, 2, 6)
db1HistEnv_Rnd = hist(db1Env_Rnd, db1Bin);
bar(db1Bin, db1HistEnv_Rnd);
xlabel('Rand TM Score Envelope'); ylabel('Count')

%Saves the figure
saveas(hFIG(iFig), fullfile(chFigDir, sprintf('%sCBASS_LFP_RandLFP.png', cBAND_LABEL{iBnd})));
%% Set the threshold for significance as what has only 5% chance of happening by chanced
% Detects significant peak: i.e. pulses
db1ThrPeak = quantile(db1Score_Rnd(bl1Peak_Rnd), .95);
bl1Pulse = bl1Peak & db1Score > db1ThrPeak;
% Detects significant epochs i.e. bouts
db1ThrEnv = quantile(db1Env_Rnd, .95);
bl1Bout = db1Env > db1ThrEnv;

iFig = iFig + 1; hFIG(iFig) = figure('Position', [50 50 1250 600]);
% Plots the histogram and the example trace of the pulses
subplot(3, 2, 1)
bar(db1Bin, [db1HistBout; db1HistCBASS_Rnd]'); hold on
db1YL = ylim; plot([1 1] * db1ThrPeak, db1YL, '--r');
legend({'LFP', 'RandData', '.05% Threshold'}, 'Location', 'NorthWest');
[~, dbP_KS] = kstest2(db1Score(bl1Peak), db1Score_Rnd(bl1Peak_Rnd));
if dbP_KS > 0.0001, chP = sprintf('%.4f', dbP_KS); else, chP = '<0.0001'; end
title(sprintf('Peaks: LFP vs Surrogate Data (KS Test: p = %s)', chP));

subplot(3, 2, 2)
db1Scr_Sel  = db1Score(in1T_Idx);
plot(db1Time, db1Scr_Sel), hold on
plot(db1Time(bl1Pulse(in1T_Idx)), db1Scr_Sel(bl1Pulse(in1T_Idx)), 'xr')
db1XL = xlim;
plot(db1XL, [1 1] * db1ThrPeak, 'r--');
title('Significant pulses');

% Plots the histogram and the the example trace of the bouts
subplot(3, 2, 3)
bar(db1Bin, [db1HistEnv; db1HistEnv_Rnd]'); hold on
db1YL = ylim; plot([1 1] * db1ThrEnv, db1YL, '--r');
legend({'LFP', 'RandData', '.05% Threshold'}, 'Location', 'NorthWest');
[~, dbP_KS] = kstest2(db1Env, db1Env_Rnd);
if dbP_KS > 0.0001, chP = sprintf('%.4f', dbP_KS); else, chP = '<0.0001'; end
title(sprintf('Envelope: LFP vs Surrogate Data (KS Test: p = %s)', chP));

subplot(3, 2, 4)
db1Scr_Sel  = db1Score(in1T_Idx);
db1Env_Sel  = db1Env(in1T_Idx);
bl1Sig_Sel  = bl1Bout(in1T_Idx);
[in1On, in1Off] = CBASS_U_FindONnOFFPoints(bl1Sig_Sel);
plot(db1Time, [db1Scr_Sel; db1Env_Sel]), hold on
for iBt = 1:length(in1On)
   plot(db1Time(in1On(iBt):in1Off(iBt)), db1Env_Sel(in1On(iBt):in1Off(iBt)), 'g')  
end
db1XL = xlim;
plot(db1XL, [1 1] * db1ThrEnv, 'r--');
title('Significant bouts');

%Creates event triggered averages of the LFP and or the CSD around bouts
    % Sets some parameters
inInterp    = 20;
in1ItrIdx   = interp1(2:inNChan - 1, 2:inNChan - 1, 2:1/inInterp:inNChan - 1);
dbScale     = .5;

    % Creates and plots the ETA
dbETA_WinSec    = [-dbCycLen dbCycLen] * 3;
in1ETA_RelIdx   = round(inSampleRate * dbETA_WinSec(1)): round(inSampleRate .* dbETA_WinSec(2));
in1EventIdx     = find(bl1Pulse);
bl1Rem          = in1EventIdx <= -in1ETA_RelIdx(1) | in1EventIdx >= inNSample - in1ETA_RelIdx(end);
in1EventIdx(bl1Rem) = [];
inNEvt          = length(in1EventIdx);
db3LFP_ETA = nan(inNChan, length(in1ETA_RelIdx), inNEvt);
for iEvt = 1:inNEvt
    db3LFP_ETA(:, :, iEvt) = db2LFP(:, in1EventIdx(iEvt) + in1ETA_RelIdx);
end
db2LFP_ETA = nanmean(db3LFP_ETA, 3);
subplot(3, 2, 5); hold on
for iChan = 1:inNChan
    plot(in1ETA_RelIdx ./ inSampleRate, (dbScale * (inNChan - iChan + 1)) + db2LFP_ETA(iChan, :), 'k');
end
title('LFP: Bout Triggered Average '); xlabel('Time (s)');

    % Computes the CSD of the ETA
db2CSD_ETA = diff(db2LFP_ETA, 2, 1);
db2CSD_ETA_Inter = interp1(2:inNChan - 1, db2CSD_ETA, in1ItrIdx);
subplot(3, 2, 6)
imagesc(in1ETA_RelIdx ./ inSampleRate, in1ItrIdx, db2CSD_ETA_Inter);
title('CSD: Bout Triggered Average'); xlabel('Time (s)'); ylabel('Channel');

saveas(hFIG(iFig), fullfile(chFigDir, sprintf('%sCBASS_Significane.png', cBAND_LABEL{iBnd})));
%% Estimates how much the identified bouts overlap with the enriched region
inFactor    = inSampleRate * 0.005;
in1EventIdx_Round = round(in1TroughIdx(in1BrchIdx == 1)'./inFactor); % Clarifies event idx
in1BoutIdx_Round  = round(find(bl1Pulse)./inFactor);

dbNewEvent  = 1 - mean(ismember(in1BoutIdx_Round, in1EventIdx_Round));
dbKeptEvent = mean(ismember(in1EventIdx_Round, in1BoutIdx_Round));

fprintf('The bout selection procedure retains %.2f%% of the original events\r', dbKeptEvent * 100);
fprintf('%.2f%% of the events are newly identified\r', dbNewEvent * 100);
%% Computes the spectra of the test
    % Sets some parameters
dbWinLenSec     = .5; % the length of the window determines the degree of smoothing
inTopFreq       = 120;

    % Computes the ETA in the window of interest
dbETA_WinSec    = [-dbWinLenSec dbWinLenSec]./2;
in1ETA_RelIdx   = round(inSampleRate * dbETA_WinSec(1)): round(inSampleRate .* dbETA_WinSec(2));
in1EventIdx     = find(bl1Pulse);
bl1Rem          = in1EventIdx <= -in1ETA_RelIdx(1) | in1EventIdx >= inNSample - in1ETA_RelIdx(end);
in1EventIdx(bl1Rem) = [];
inNEvt          = length(in1EventIdx);
db3LFP_ETA = nan(inNChan, length(in1ETA_RelIdx), inNEvt);
for iEvt = 1:inNEvt
    db3LFP_ETA(:, :, iEvt) = db2LFP(:, in1EventIdx(iEvt) + in1ETA_RelIdx);
end
db2LFP_ETA = nanmean(db3LFP_ETA, 3);

% Computes the spectrum of the ETA and plots it
db2FFT      = fft(db2LFP_ETA .* hamming(size(db2LFP_ETA, 2))', [], 2);
db2FFT      = db2FFT(:, 1:length(in1ETA_RelIdx)/2+1);
db2Power    = (1/dbWinLenSec).*abs(db2FFT).^2;
db2Power(:, 2:end-1) = 2*db2Power(:, 2:end-1);

% Makes the figure
iFig = iFig + 1; hFIG(iFig) = figure;
db1Freq = 0:1/dbWinLenSec:inTopFreq;
plot(db1Freq, 10 * log10(mean(db2Power(:, 1:length(db1Freq)))));
%% Attempt to see what happens when the animal is running
in1EventIdx = in1TroughIdx(in1BrchIdx == 1);
bl1Event    = false(size(bl1Pulse)); bl1Event(in1EventIdx) = true;

% Calculates the frequency of bouts and events during quiet and running 
dbFreqBoutRun   = inSampleRate * sum(bl1Pulse & bl1Run)./sum(bl1Run);
dbFreqBoutQuiet = inSampleRate * sum(bl1Pulse & ~bl1Run)./sum(~bl1Run);
dbFreqEvntRun   = inSampleRate * sum(bl1Event & bl1Run)./sum(bl1Run);
dbFreqEvntQuiet = inSampleRate * sum(bl1Event & ~bl1Run)./sum(~bl1Run);
% Relative increase of frequency
dbRatioBout = dbFreqBoutRun./dbFreqBoutQuiet;
dbRatioEvnt = dbFreqEvntRun./dbFreqEvntQuiet;
% Significance of the increase
dbPVal_Bout =  2 * (1 - normcdf(abs(dbFreqBoutRun - dbFreqBoutQuiet)./...
    sqrt((dbFreqBoutRun ./ sum(bl1Pulse & bl1Run)) +  (dbFreqBoutQuiet ./ sum(bl1Pulse & ~bl1Run)) )));
dbPVal_Evnt =  2 * (1 - normcdf(abs(dbFreqEvntRun - dbFreqEvntQuiet)./...
    sqrt((dbFreqEvntRun ./ sum(bl1Event & bl1Run)) +  (dbFreqEvntQuiet ./ sum(bl1Event & ~bl1Run)) )));

% Print the results on the screen
if dbPVal_Evnt > 0.0001, chP = sprintf('%.4f', dbPVal_Evnt); else, chP = '<0.0001'; end
fprintf('Original Troughs:\tQuiet: %.2fHz\tRun: %.2fHz\tRatio:%.2f\t(p = %s)\r', ....
    dbFreqEvntQuiet, dbFreqEvntRun, dbRatioEvnt, chP);
if dbPVal_Bout > 0.0001, chP = sprintf('%.4f', dbPVal_Bout); else, chP = '<0.0001'; end
fprintf('Selected Bouts:\t\tQuiet: %.2fHz\tRun: %.2fHz\tRatio:%.2f\t(p = %s)\r', ....
    dbFreqBoutQuiet, dbFreqBoutRun, dbRatioBout, chP);

% Plots moving averages of the bouts and the initial events
in1RunON = 1 + find(~bl1Run(1:end-1) & bl1Run(2:end));

% Convolutes the trace by a square window to compute an instantaneous
% frequency
dbWinLenSec     = 0.5;
db1BoutTrace    = conv(bl1Pulse, rectwin(dbWinLenSec * inSampleRate), 'same')./dbWinLenSec;
db1EventTrace   = conv(bl1Event, rectwin(dbWinLenSec * inSampleRate), 'same')./dbWinLenSec;

% Defines the ETA 
dbETA_WinSec    = [-5 10];
in1ETA_RelIdx   = round(inSampleRate * dbETA_WinSec(1)): round(inSampleRate .* dbETA_WinSec(2));
db1ETA_Bout     = mean(db1BoutTrace(in1ETA_RelIdx + in1RunON'));
db1ETSEM_Bout   = std(db1BoutTrace(in1ETA_RelIdx + in1RunON'))./sqrt(length(in1RunON));
db1ETA_Event    = mean(db1EventTrace(in1ETA_RelIdx + in1RunON'));
db1ETSEM_Event  = std(db1EventTrace(in1ETA_RelIdx + in1RunON'))./sqrt(length(in1RunON));

% Defines the time
db1Time     = in1ETA_RelIdx./inSampleRate;

% Plot the results
iFig = iFig + 1; hFIG(iFig) = figure('Position', [50 50 1250 600]);
subplot(1, 2, 1)
fill([db1Time db1Time(end:-1:1)], ...
    [db1ETA_Event + db1ETSEM_Event  db1ETA_Event(end:-1:1) - db1ETSEM_Event(end:-1:1)], ...
    [.5 .5 .5], 'LineStyle', 'none', 'FaceAlpha', .3); hold on
plot(db1Time, db1ETA_Event, 'k'); hold on
db1YL = ylim; plot([0 0], db1YL, 'r--');
ylabel('Instantaneous Frequency (Hz)'); xlabel('Time (s)');
title('Original Event Around Running Onset');
subplot(1, 2, 2)
fill([db1Time db1Time(end:-1:1)], ...
    [db1ETA_Bout + db1ETSEM_Bout  db1ETA_Bout(end:-1:1) - db1ETSEM_Bout(end:-1:1)], ...
    [.5 .5 .5], 'LineStyle', 'none', 'FaceAlpha', .3); hold on
plot(db1Time, db1ETA_Bout, 'k'); hold on
% plot(db1Time, (db1ETA_Bout + [db1ETSEM_Bout; -db1ETSEM_Bout])', 'k--');
db1YL = ylim; plot([0 0], db1YL, 'r--');
ylabel('Instantaneous Frequency (Hz)'); xlabel('Time (s)');
title('Selected Bouts Around Running Onset');

saveas(hFIG(iFig), fullfile(chFigDir, sprintf('%s_Pulse_Running.png', cBAND_LABEL{iBnd})));
%% Plots exemple of the traces around bout onset
dbEx_WinSec     = [-2 3];
in1Ex_RelIdx    = round(inSampleRate * dbEx_WinSec(1)): round(inSampleRate .* dbEx_WinSec(2));
db1Time         = in1Ex_RelIdx ./ inSampleRate;
dbScale = 2;

for iPlt = 1:9
    iFig = iFig + 1; hFIG(iFig) = figure('Position', [50 50 1250 600]); hold on
    in1Sel = in1Ex_RelIdx + in1RunON(iPlt);    
    for iChan = 1:inNChan
        plot(db1Time, (dbScale * (inNChan - iChan + 1)) + db2LFP(iChan, in1Sel), 'k');
    end
    plot(db1Time, db1Score(in1Sel)); 
    in1Pulse = find(bl1Pulse(in1Sel));
    db1YL = ylim;
    for iPls = 1:length(in1Pulse)
        plot([1 1] * db1Time(in1Pulse(iPls)), db1YL, '--r')
    end
    plot([0 0], db1YL, 'g', 'LineWidth', 2)
    title('Example Pulses (red) around running onset (green)'); 
    xlabel('Time (s)');
    saveas(hFIG(iFig), fullfile(chFigDir, sprintf('%s_Pulse_Example_%d_Running.png', cBAND_LABEL{iBnd}, iPlt)));
end