% Hi my friend here is a little script to 1. Load some exemple data that
% you can play with  and 2. show you an exemple of what I have been pocking
% around with with generalized eigen value decomposition to identify filters. 
% I hope that you will find it interesting.

%Set the current directory to where the script is
clear, close all
% chMyAdress = '/media/storage/Quentin/Collaboration/Antonio/2020-02-04/'; % < -------------Edit 
chMyAdress = '/gpfs/ysm/scratch60/dietrich/ahf38/Cardin/quentin/'; % < -------------Edit 
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
db2LFP = sINPUT_2.sL4MMDS.db2LFP; %Exctracts the LFP matrix
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
%% Example of an attempt to do some generalized eigen value decomposition of the spatio temporal propagation of activity
%Focuses on the presentation and filters the LFP in the range of interest
% db1Band = [15 30]; % Beta (i.e. the rythm induced by visual stimulation)
db1Band = [30 80]; % Gamma (i.e the rythm induced by running ... and licking but there is not licking data here)
[B, A] = butter(2, 2*db1Band/inSampleRate); 

% Filters the LFP
db2_Filt_LFP = filtfilt(B, A, db2LFP')';

%Defines widow length and defines the size of overlapping chunks of the
%data. This allows to subsample the signal while still detectin a potential
%propagation of activity. I am not totally sure that it is the best way but
%it is a first pass
dbWinLenSec = 3/db1Band(1); %About five cycles of the lower edge of the frequency of interest
inNStep     = 1; % Adjacent windows with no overlap. 2 for 2 window onset per window length. 4 for 4 windows onset per window length (must divide the window length evenly)
inWinLen    = round(dbWinLenSec * inSampleRate);
inStepLen   = inWinLen / inNStep;

%Creates a reference trace
db1Trace = db2_Filt_LFP(1, :);

%Checks for the maximal number of chunks that can fit the trace and removes
%the points that might not be used
inNChunk = floor(length(db1Trace) / inStepLen);
% db1Trace = db1Trace(1:inNChunk * inStepLen); %Removes unused points
inNChunk = inNChunk - (inNStep - 1);

%Create a matrix, whose columns are samples of the signal whose length is
%windowLengthSec. Successive samples (columns) are overlapping and spaced
%by windowLength/step seconds. 
in1ChunkBegIdx  = ((1:inNChunk) * inStepLen) - (inStepLen - 1);
in1ChunkIdx     = ((1:inWinLen) - 1)';
in2ChunkIdx     = repmat(in1ChunkBegIdx, inWinLen, 1) + repmat(in1ChunkIdx, 1, inNChunk);

%Gets the central index of each window
in1ChunkIdx = ceil(median(in2ChunkIdx));

%Computes the indices of exclusive conditions
bl1RunOnly_Chnk      = bl1RunOnly(in1ChunkIdx);
bl1PresOnly_Chnk     = bl1PresOnly(in1ChunkIdx);
bl1Quiet_Chnk        = bl1Quiet(in1ChunkIdx);

%Construct a chunk map of the LFP. 
db2LFP_ST_Chunks = [];
%in1ChanSel = [3 6 9 12]; % Selects a subset of the channels to save time
in1ChanSel = 1:size(db2_Filt_LFP, 1);
inNChan = length(in1ChanSel);
for iChan = in1ChanSel
    db1Trace = db2_Filt_LFP(iChan, :);
%     db2LFP_ST_Chunks = db1Trace(in2ChunkIdx);
    db2LFP_ST_Chunks = cat(1, db2LFP_ST_Chunks , db1Trace(in2ChunkIdx));
end

%Computes the spatio temporal covariance matrix of the activity for all
%conditions
db2COV_Run      = cov(db2LFP_ST_Chunks(:, bl1RunOnly_Chnk)');
db2COV_Pres     = cov(db2LFP_ST_Chunks(:, bl1PresOnly_Chnk)');
db2COV_Quiet    = cov(db2LFP_ST_Chunks(:, bl1Quiet_Chnk)');
% figure, imagesc(db2COV_Run); title('Run');
% figure, imagesc(db2COV_Pres); title('Pres');
% figure, imagesc(db2COV_Quiet); title('Quiet');
fprintf('Done\r')

%Computes the generalized eigenvalue decomposition for running and
%presentation
% tic, [db2EigVec_Pres, db2EigVal_Pres]    = eig(db2COV_Pres, db2COV_Quiet); toc
% % takes the real part (assumes that complexity arises because of round offs errors which would need to be verified -- interpret with caution)
% db2EigVec_Pres = real(db2EigVec_Pres); 
% db2EigVal_Pres = real(db2EigVal_Pres);
tic, [db2EigVec_Run, db2EigVal_Run]      = eig(db2COV_Run, db2COV_Quiet); toc
% takes the real part (assumes that complexity arises because of round offs errors which would need to be verified -- interpret with caution)
db2EigVec_Run = real(db2EigVec_Run);
db2EigVal_Run = real(db2EigVal_Run);
%% Reconstitute a forward model to compute filters and plots the 10 best filters (strongest eigenvalues)
inNCmp = 10;

% figure('Position', [100 100 1600 900])
% for iCmp = 1:inNCmp
%     subplot(1, inNCmp, iCmp)
%     db1Fwrd = db2COV_Pres * db2EigVec_Pres(:, end - iCmp + 1);
%     db2Fwrd = reshape(db1Fwrd, inWinLen, inNChan)';
%     imagesc(db2Fwrd), title(sprintf('Pres Cmp %d', iCmp))
% ends

figure('Position', [100 100 1600 900])
for iCmp = 1:inNCmp
    subplot(1, inNCmp, iCmp)
    db1Fwrd = db2COV_Run * db2EigVec_Run(:, end - iCmp + 1);
    db2Fwrd = reshape(db1Fwrd, inWinLen, inNChan)';
    imagesc(db2Fwrd), title(sprintf('Run Cmp %d', iCmp))
end

% Commments: So here I find the result somewhat interesting but the
% selection of filters is where I need to explore. I hope you will find it
% interesting. Let's talk soon.
