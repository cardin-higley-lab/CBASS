% Hi my friend here is a little script to 1. Load some exemple data that
% you can play with  and 2. show you an exemple of what I have been pocking
% around with with generalized eigen value decomposition to identify filters. 
% I hope that you will find it interesting.

%Set the current directory to where the script is
clear, close all
% addpath('/gpfs/ysm/home/ahf38/Documents/gamma_bouts/scripts')
% chMyAdress = '/media/storage/Quentin/Collaboration/Antonio/2020-02-04/'; % < -------------Edit 
chMyAdress = 'D:\gamma_bouts\scripts'; % < -------------Edit 
% chMyAdress = '/gpfs/ysm/home /ahf38/Documents/gamma_bouts/data/Example_1'; % < -------------Edit 
cd(chMyAdress);
%% Load and formats the exemple recording
%Recreates the session folder name
chSessionName       = 'Example_1';

%Checks the for visual stimulation meta structure
chDPSFile = 'PCP_L1_DetectPresentationSet.mat';
%Checks the for ephys meta structure
chMDSFile = 'PCP_L4_MakeMetaDataStructure.mat';

%Loads the input
sINPUT_1 = load(fullfile(chSessionName, chDPSFile), '-mat'); sINPUT_1 = sINPUT_1.sCFG;
sINPUT_2 = load(fullfile(chSessionName, chMDSFile), '-mat'); sINPUT_2 = sINPUT_2.sCFG;

%Get the overal sample rate
inSampleRate = sINPUT_2.sL4MMDS.inWorkSampleRate;

%Extracts the LFP
db2LFP = sINPUT_2.sL4MMDS.db2LFP; %Exctracts the LFP matrixgit
db2LFP = PCP_U_NormalizeLFP(db2LFP, inSampleRate); % Filters and z-scores the LFP 
db2LFP = db2LFP(2:end, :); % Removes the first channel of the LFP which is set as the reference. We are left with 15 channels.

%Gets the indices of when a visual presentation was shown
in1PresOnIdx    = CBASS_U_GetTStampEventIndex(sINPUT_2.sL4MMDS.db1TStamps, sINPUT_1.sL1DP.db1PresOnTStamp);
in1PresOffIdx   = CBASS_U_GetTStampEventIndex(sINPUT_2.sL4MMDS.db1TStamps, sINPUT_1.sL1DP.db1PresOffTStamp);
bl1Pres         = CBASS_U_MakeEpochVector(in1PresOnIdx, in1PresOffIdx, length(sINPUT_2.sL4MMDS.db1TStamps));

%Gets the indices of when the mouse is running or whisking
bl1Run          = sINPUT_2.sL4MMDS.bl1WheelOn;

%Computes the indices of exclusive conditions
bl1RunOnly      = bl1Run & ~bl1Pres;
bl1PresOnly     = ~bl1Run & bl1Pres;
bl1Quiet        = ~bl1Run & ~bl1Pres;
%% performs PCA before hand to demix the LFP the phase randomization then recompose the signal
[db2Coeff, db2LFP_Rnd] = pca(db2LFP');
db2LFP_Rnd = CBASS_U_PhaseRandomize1D(db2LFP_Rnd, 1);
db2LFP_Rnd = (db2LFP_Rnd/db2Coeff)';
%% Compares the spectra of the original LFP and of the phase randomized LFP
dbWinLenSec = .5; %Default of CBASS_U_FourierPower
inTopFreq   = 120; %Maximum frequency to plot in Hz
in1FrqIdx   = 1:inTopFreq*dbWinLenSec+1;
in1FrqLbl   = 0:1/dbWinLenSec:inTopFreq;

inNChan = size(db2LFP, 1);
dbScale = 3;
    
[db2Pwr, db2PwrRnd] = deal([]);
for iChan = 1:inNChan
    % Aggregates the power for the LFP
    db2Power = CBASS_U_FourierPower(db2LFP(iChan, :), inSampleRate);
    db2Pwr = cat(2, db2Pwr, mean(db2Power, 2));
    
    % Aggregates the power for the phase randomized LFP
    db2Power = CBASS_U_FourierPower(db2LFP_Rnd(iChan, :), inSampleRate);
    db2PwrRnd = cat(2, db2PwrRnd, mean(db2Power, 2));
end

figure
% Plots the spectra
hPLT(1) = subplot(2, 2, 1); plot(in1FrqLbl, 10 * log10(mean(db2Pwr(in1FrqIdx, :), 2)));
title('LFP Spectrum'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
hPLT(2) = subplot(2, 2, 2); plot(in1FrqLbl, 10 * log10(mean(db2PwrRnd(in1FrqIdx, :), 2)));
title('Phase Randomized LFP Spectrum'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
in1PltIdx   = 1:10*inSampleRate; 
in1PltLabel = in1PltIdx / inSampleRate;
linkaxes(hPLT);
% Plots the 10 first seconds of the LFP and the phase randomized LFP
hPLT2(1) = subplot(2, 2, 3); hold on
for iChan = 1:inNChan
    plot(in1PltLabel, (dbScale * (inNChan - iChan + 1)) + db2LFP(iChan, in1PltIdx), 'k');
end
title('LFP Exerpt'); xlabel('Time (s)'); 
hPLT2(2) = subplot(2, 2, 4); hold on
for iChan = 1:inNChan
    plot(in1PltLabel, (dbScale * (inNChan - iChan + 1)) + db2LFP_Rnd(iChan, in1PltIdx), 'k');
end
title('Phase Randomized LFP Exerpt'); xlabel('Time (s)');
linkaxes(hPLT2);
%% Formats the LFP to extract beta and gamma propagation motifs
%Defines Beta (15-30Hz i.e. the rythm induced by visual stimulation) and
%Gamma % Gamma (30-80Hz i.e the rythm induced by running ... and licking
%but there is not licking data here)
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND = {[15 30], [30 80]};

%Initalizes the motifs cell array
[cMOTIF, cTROUGH_IDX]  = deal(cell(size(cBAND)));

%Defines a reference channel (Here we are taking a channel that is roughly
%close to layer 4.
inRefChan = 5;

%Loops over bands
for iBnd = 1:length(cBAND)
    % Filters the LFP
    [B, A] = butter(2, 2 * cBAND{iBnd} / inSampleRate); 
    db2_Filt_LFP = filtfilt(B, A, db2LFP_Rnd'); %The LFP is transposed because filtfilt and hilbert opperates over rows

    % Uses the Hilbert Transform
    db2_Hilbert = hilbert(db2_Filt_LFP);
    
    % Computes the amplitude and the phase
    db1_Amp     = abs(db2_Hilbert); 
    db1_Phase   = angle(db2_Hilbert);
    
    % Finds the indices of troughs on the reference channels
    bl1RefTrough = db1_Phase(1:end-1, inRefChan) < 0 & db1_Phase(2:end, inRefChan) > 0;
    
    % Formats the bouts so that that each row is a motif (the ref channel
    % is at the trough of the oscillation) and that col 1:15 correspond to
    % the amplitude of each channel and col 16:30 correspond to the phase.
    % Note these two groups should be normalized separately
    cMOTIF{iBnd}        = [db1_Amp(bl1RefTrough, :) db1_Phase(bl1RefTrough, :)];
    cTROUGH_IDX{iBnd}   = find(bl1RefTrough);
end

%% Tries to see if preferential regions of the motifs space coincide with specific states
for iBnd = 1:length(cBAND)
    % Computes a pca of the motifs    
    [~, db1PC] = pca(cMOTIF{iBnd});
    
    % Plots the PCA with labels
    figure('Position', [50 50 1250 600])
    
    % Visual Stimulation
    subplot(1, 2, 1);
    gscatter(db1PC(:,1), db1PC(:,2), bl1Pres(cTROUGH_IDX{iBnd}));
    xlabel('PC1'); ylabel('PC2');
    title(sprintf('%s (0 = No Pres; 1 = Pres)', cBAND_LABEL{iBnd}))
    
    % Running
    subplot(1, 2, 2);
    gscatter(db1PC(:,1), db1PC(:,2), bl1Run(cTROUGH_IDX{iBnd}));
    xlabel('PC1'); ylabel('PC2');
    title(sprintf('%s (0 = No Running; 1 = Running)', cBAND_LABEL{iBnd}))

end