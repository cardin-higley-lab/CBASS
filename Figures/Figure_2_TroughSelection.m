%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
iExp = 2;
chExperiment = cEXP{iExp};
chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/'; 
chDir = fullfile(chPath,'Data', chExperiment);
% chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath, 'Figures', 'Figure_2_TroughSelection');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Loads the data
sREC    = CBASS_L0_LoadData(chDir);
% sREC    = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Sets general parameter
% Sets the band and the state of interest
chLabel         = 'Beta';
db1Band         = [15 30];
bl1State        = sREC.bl1Pres;
chStateLabel    = 'Stim';

% Sets options for trough extraction
chDataFormat    = 'complex';
blZScore        = true;
inRefChan       = 5;

% Sets parameters for plotting the trough
in1StimON   = find(~bl1State(1:end - 1) & bl1State(2:end));
inAnchor    = in1StimON(2) - (0.2 * sREC.inSampleRate);
db1WinSec   = [0 .6];

% Sets the gloabal label for the experiment;
chExpLbl    = [chExperiment '_' chLabel];

% Gets the troughs
sTROUGH = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);

%% Initializizes the figure
clear hFIG
hFIG = figure('Position', [50 50 1250 600]);
cFIG_NAME = {chExpLbl};

% Plots the pulse triggerred averages
subplot(1, 3, 1);
CBASS_Plot_LFP_EventExample(sREC.db2LFP, inAnchor, db1WinSec, sREC.inSampleRate);
title('LFP');

% Plots event enrichment
% Filters the LFP
[B, A] = butter(2, 2 * db1Band / sREC.inSampleRate);
db2_Filt_LFP = filtfilt(B, A, sREC.db2LFP')'; %The LFP is transposed because filtfilt and hilbert opperates over rows
subplot(1, 3, 2);
CBASS_Plot_LFP_EventExample(db2_Filt_LFP, inAnchor, db1WinSec, sREC.inSampleRate, {sTROUGH.in1Index}, {'Troughs'});
title(sprintf('%d-%dHz Filtered LFP', db1Band(1), db1Band(2)))
% Appends the reference channel
    % computes the scaling
db1YL       = ylim; db1XL       = xlim;
inNChan     = size(db2_Filt_LFP, 1);
dbScale     = db1YL(2) ./ (inNChan + 1.5);
    % Selects the wanted chunk of the LFP
in1RelIdx   = round(sREC.inSampleRate * db1WinSec(1)): round(sREC.inSampleRate .* db1WinSec(2));
db1Time     = in1RelIdx ./ sREC.inSampleRate;
in1Sel      = in1RelIdx + inAnchor;
% Plots the reference channel
 plot(db1Time, (dbScale * (inNChan - inRefChan + 1)) + db2_Filt_LFP(inRefChan, in1Sel), 'r', 'DisplayName', 'Ref');
 
% Plots the inset around a trough
    % Calculates the coordinates the inset
in1T_Anchor = sTROUGH.in1Index(ismember(sTROUGH.in1Index, in1Sel));
inT_Anc     = in1T_Anchor(5);
in1Idx_Bnd  = round(sREC.inSampleRate .* [-.5 .5] ./ median(db1Band));
in1RelIdx_2 = in1Idx_Bnd(1):in1Idx_Bnd(2);
db1Time_2   = in1RelIdx_2 ./ sREC.inSampleRate;
in1Sel_2    = in1RelIdx_2 + inT_Anc;
db1TrghEx   = db2_Filt_LFP(inRefChan, in1Sel_2);
db1YL       = [min(db1TrghEx) - (range(db1TrghEx) .* .3) max(db1TrghEx) + (range(db1TrghEx) .* .3)];
    % Calculates the coordinates of the square in the third plot
db1SqrX     = (in1Sel_2([1 end end 1 1]) - in1Sel(1) - 1) ./ sREC.inSampleRate;
db1SqrY     = db1YL([1 1 2 2 1]) + dbScale * (inNChan - inRefChan + 1);
plot(db1SqrX, db1SqrY, '--r', 'LineWidth' , 1, 'DisplayName', 'Inset')

subplot(1, 3, 3)
plot(db1Time_2, db1TrghEx, 'r'); hold on
plot([0 0], db1YL, '--k');
title('Example Trough'); ylim(db1YL), xlim(db1Time_2([1 end]));
xlabel('Time (s)'); ylabel('Filtered LFP (mV)')

%% Saves the figure
CBASS_SaveFig(chOutPath, hFIG, cFIG_NAME);
close all;