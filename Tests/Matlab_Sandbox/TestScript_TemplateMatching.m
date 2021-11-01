% Runs the data analysis script
'D:\gamma_bouts\scripts'; % < -------------Edit
run LoadDataScript_BetaGamma_HilbertFormatting.m
close all
load('ClusterBranchIndices_v2.mat');

% Creates the figure directory
chFigDir = 'TemplateMatching';
if ~exist(chFigDir, 'dir'); mkdir(chFigDir); end
iFig = 0;
%%
% Gets the indices of the epoches on the initial LFP recording
in1TroughIdx = cTROUGH_IDX{2};

% Gets the size of the LFP
[inNChan, inNSample] = size(db2LFP);

% Gets the indices of the branches
in1BrchIdx  = in1BranchIdx_40Clusters;
inNBch      = max(in1BrchIdx);

% Sets the indices of the event triggered average
dbETA_WinSec    = [-.05 .05];
in1ETA_RelIdx   = round(inSampleRate * dbETA_WinSec(1)): round(inSampleRate .* dbETA_WinSec(2));

% Initalizes a figure and sets interpolation for the CSD and the scale for
% plotting LFP averages
iFig = iFig + 1; hFIG(iFig) = figure('Position', [50, 50, 1250, 600]);
inInterp    = 20;
in1ItrIdx   = interp1(2:inNChan - 1, 2:inNChan - 1, 2:1/inInterp:inNChan - 1);
dbScale     = .5;

for iBch = 1:inNBch
    % Finds the trigger events (the LFP band trouch of the selected branch)
    % and removes them if they too close to the edge of the recording
    in1EventIdx     = in1TroughIdx(in1BrchIdx == iBch);
    bl1Rem          = in1EventIdx <= -in1ETA_RelIdx(1) | in1EventIdx >= inNSample - in1ETA_RelIdx(end);
    in1EventIdx(bl1Rem) = [];
    inNEvt          = length(in1EventIdx);
    
    % Creates the ETA
    db3LFP_ETA = nan(inNChan, length(in1ETA_RelIdx), inNEvt);
    for iEvt = 1:inNEvt
        db3LFP_ETA(:, :, iEvt) = db2LFP(:, in1EventIdx(iEvt) + in1ETA_RelIdx);
    end
    db2LFP_ETA = nanmean(db3LFP_ETA, 3);
    
    % Computes the CSD of the ETA
    db2CSD_ETA = diff(db2LFP_ETA, 2, 1);
    db2CSD_ETA_Inter = interp1(2:inNChan - 1, db2CSD_ETA, in1ItrIdx);
    
    % Plots the figure
    % LFP
    subplot(2, inNBch, iBch); hold on
    for iChan = 1:inNChan
        plot(in1ETA_RelIdx * inSampleRate, (dbScale * (inNChan - iChan + 1)) + db2LFP_ETA(iChan, :), 'k');
    end
    title(sprintf('LFP: 40 Clusters: Branch %d ', iBch)); xlabel('Time (s)');
    % CSD
    subplot(2, inNBch, inNBch + iBch)
    imagesc(in1ETA_RelIdx, in1ItrIdx, db2CSD_ETA_Inter);
    title('CSD ETA'); xlabel('Time (s)'); ylabel('Channel');
end

saveas(hFIG(iFig), fullfile(chFigDir, 'ETA_LFP_CSD_40Clusters.png'))
%% Performs a match algorithm on the trace
% Cleans workspace
clear sOPTION

% Defines Bands
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND = {[15 30], [30 80]};

% Defines what band we are dealing with here
iBnd = 2;

% Gets the size of the LFP
[inNChan, inNSample] = size(db2LFP);

% Gets the indices of the branches
in1BrchIdx  = in1BranchIdx_40Clusters;

% Gets the indices of the epoches on the initial LFP recording
in1TroughIdx    = cTROUGH_IDX{2};
in1EventIdx     = in1TroughIdx(in1BrchIdx == 1); % Clarifies event idx
inNEvt          = length(in1EventIdx);

% Sets the normalization options
cCNTR_NOMR      = {'Raw', 'Centered', 'Normalized', 'CenterNormalized'};
bl1Center       = [false true false true];
bl1Normalize    = [false false true true];

% Loops over normalization condition
for iCnNr = 1:length(bl1Center)
    % Sets whether the scores should be centered about the mean and normalized
    % by the moving average of the signal
    blCntr  = bl1Center(iCnNr);
    blNrm   = bl1Normalize(iCnNr);
    
    % Option 1 ----------------------------------------------------------------
    % Takes the dot product of the Analytic signal yielded by the Hilbert
    % transform and the average events (i.e. the trough of enriched clusters)
    % in the complex plane
    
    % Filters the LFP
    [B, A] = butter(2, 2 * cBAND{iBnd} / inSampleRate);
    db2_Filt_LFP = filtfilt(B, A, db2LFP'); %The LFP is transposed because filtfilt and hilbert opperates over rows
    
    % Uses the Hilbert Transform
    db2_Hilbert = hilbert(db2_Filt_LFP);
    
    % Defines the filter
    db1Filter = mean(db2_Hilbert(in1TroughIdx, :));
    
    % Stores the output
    sOPTION(iCnNr, 1).db1RawScore  = CBASS_U_MultiChannelTemplateMatching(db2_Hilbert.', db1Filter.', blCntr, blNrm);
    sOPTION(iCnNr, 1).db1Envelope  = abs(hilbert(sOPTION(iCnNr, 1).db1RawScore));
    sOPTION(iCnNr, 1).chName       = 'ComplexTrough';
    
    % Option 2 ----------------------------------------------------------------
    % Separates the amplitude and phase of the analytical signal and does the
    % templated matching there.
    db2Signal = [abs(db2_Hilbert) angle(db2_Hilbert)];
    % Defines the filter
    db1Filter = mean(db2Signal(in1TroughIdx, :));
    
    % Stores the output
    sOPTION(iCnNr, 2).db1RawScore  = CBASS_U_MultiChannelTemplateMatching(db2Signal', db1Filter', blCntr, blNrm);
    sOPTION(iCnNr, 2).db1Envelope  = abs(hilbert(sOPTION(iCnNr, 2).db1RawScore));
    sOPTION(iCnNr, 2).chName       = 'PolarTrough';
    
    % Option 3 ----------------------------------------------------------------
    % Crosscorrelation of the raw LFP and the event triggered average 1 cycle
    
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
    
    % Stores the output
    sOPTION(iCnNr, 3).db1RawScore  = CBASS_U_MultiChannelTemplateMatching(db2LFP, db2Filter, blCntr, blNrm);
    sOPTION(iCnNr, 3).db1Envelope  = abs(hilbert(sOPTION(iCnNr, 3).db1RawScore));
    sOPTION(iCnNr, 3).chName       = 'TempMatch';
    
    % Option 4 ----------------------------------------------------------------
    % Crosscorrelation of the raw LFP and the event triggered average 1 cycle
    sOPTION(iCnNr, 4).db1RawScore  = CBASS_U_MultiChannelTemplateMatching(diff(db2LFP, 2, 1), diff(db2Filter, 2, 1), blCntr, blNrm);
    sOPTION(iCnNr, 4).db1Envelope  = abs(hilbert(sOPTION(iCnNr, 4).db1RawScore));
    sOPTION(iCnNr, 4).chName       = 'CSD-TM';
end
%% Calculates and plots the correlation between each option and the triggers
in1Event    = false(1, size(db2LFP, 2));
in1Event(in1EventIdx) = true;

db2CorrIm     = zeros(size(sOPTION));
for iCnNr = 1:size(sOPTION, 1)
    for iOpt = 1:size(sOPTION, 2)
        db2Corr = corr([in1Event; sOPTION(iCnNr, iOpt).db1RawScore]');
        db2CorrIm(iCnNr, iOpt) = db2Corr(1, 2);
    end
end

% Plots the figure
iFig = iFig + 1; hFIG(iFig) = figure;
imagesc(1:4, 1:4, db2CorrIm);
colorbar
set(gca, 'XTick', 1:4, 'XTickLabel', {sOPTION(1,:).chName}, ...
    'YTick', 1:4, 'YTickLabel', cCNTR_NOMR);
title('Correlation with manifold events')

saveas(hFIG(iFig), fullfile(chFigDir, 'ScoreComparision.png'));
%%  Plots the exemples of interest
clear hAX hOP_PLT

% Computes time selection for plot.
dbTimeSec   = 5;
dbTOffset   = 1215;
in1T_Idx    = (1:round(inSampleRate * dbTimeSec)) + round(inSampleRate * dbTOffset);
db1Time     = in1T_Idx ./ inSampleRate;

% Computes the events to plot
in1PltEvnt  = in1EventIdx(ismember(in1EventIdx, in1T_Idx)) - round(inSampleRate * dbTOffset);

% Sets the selection of plots to
in1CnNr_Sel     = [1 4];
in1Option_Sel   = [1, 3];

for iSel = 1:length(in1CnNr_Sel)
    % Select the option to plot
    iCnNr   = in1CnNr_Sel(iSel);
    iOpt    = in1Option_Sel(iSel);
    
    iFig = iFig + 1; hFIG(iFig) = figure('Position', [50 50 1250 600]);
    hAX(1) = subplot(4, 1, 1:3); hold on
    dbScale = 2;
    for iChan = 1:inNChan
        plot(db1Time, (dbScale * (inNChan - iChan + 1)) + db2LFP(iChan, in1T_Idx), 'k');
    end
    db1YL = ylim;
    for iEvt = 1:length(in1PltEvnt); plot(db1Time(in1PltEvnt(iEvt)) * [1 1], db1YL, 'r'); end
    hAX(2) = subplot(4, 1, 4); hold on
    hOP_PLT(1) = plot(db1Time, sOPTION(iCnNr, iOpt).db1RawScore(in1T_Idx));
    hOP_PLT(2) = plot(db1Time, sOPTION(iCnNr, iOpt).db1Envelope(in1T_Idx));
    legend(hOP_PLT, {sOPTION(iCnNr, iOpt).chName 'Envelope'});
    linkaxes(hAX, 'x');
    
    chFigName = sprintf('Example_%s_%s', cCNTR_NOMR{iCnNr}, sOPTION(iCnNr, iOpt).chName);
    saveas(hFIG(iFig), fullfile(chFigDir, [chFigName '.png']));
end