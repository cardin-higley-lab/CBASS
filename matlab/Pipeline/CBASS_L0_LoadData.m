function sREC = CBASS_L0_LoadData(chDir)
% L0 of the bout pipeline. Loads the data and return them into a structure
% sREC.

% Makes sure that utilities are added to the paths
addpath(genpath(fullfile('..','Pipeline')))
% addpath('../Pipeline/Utilities')

% Check the input
    % Check that chAdress is a directory
if ~isfolder(chDir), error('%s is not a directory', chDir); end
    %Check that the file for visual stimulation meta structure exists
chDPSFile = 'PCP_L1_DetectPresentationSet.mat';
if ~exist(fullfile(chDir, chDPSFile), 'file'), fprintf('%s not found', chDPSFile); blPres = false; 
else, blPres = true; end
    %Check that the file for ephys meta structure exists
chMDSFile = 'PCP_L4_MakeMetaDataStructure.mat';
if ~exist(fullfile(chDir, chMDSFile), 'file'), error('%s not found', chMDSFile); end

%Load the input
if blPres, sINPUT_1 = load(fullfile(chDir, chDPSFile), '-mat'); sINPUT_1 = sINPUT_1.sCFG;  end
sINPUT_2 = load(fullfile(chDir, chMDSFile), '-mat'); sINPUT_2 = sINPUT_2.sCFG;

%Get the overal sample rate
inSampleRate = sINPUT_2.sL4MMDS.inWorkSampleRate;

%Extract the LFP
db2LFP = sINPUT_2.sL4MMDS.db2LFP; %Exctracts the LFP matrix
db2LFP = PCP_U_NormalizeLFP(db2LFP, inSampleRate); % Filters and z-scores the LFP 
db2LFP = db2LFP(2:end, :); % Removes the first channel of the LFP which is set as the reference. We are left with 15 channels.

%Get the indices of when a visual presentation was shown
if blPres
    in1PresOnIdx    = CBASS_U_GetTStampEventIndex(sINPUT_2.sL4MMDS.db1TStamps, sINPUT_1.sL1DP.db1PresOnTStamp);
    in1PresOffIdx   = CBASS_U_GetTStampEventIndex(sINPUT_2.sL4MMDS.db1TStamps, sINPUT_1.sL1DP.db1PresOffTStamp);
    bl1Pres         = CBASS_U_MakeEpochVector(in1PresOnIdx, in1PresOffIdx, length(sINPUT_2.sL4MMDS.db1TStamps));
else
    bl1Pres = false(1, size(db2LFP, 2));
end

%Get the indices of when the mouse is running or whisking
bl1Run          = sINPUT_2.sL4MMDS.bl1WheelOn;

% %Compute the indices of exclusive conditions
% bl1RunOnly      = bl1Run & ~bl1Pres;
% bl1PresOnly     = ~bl1Run & bl1Pres;
% bl1Quiet        = ~bl1Run & ~bl1Pres;

% Store variables in the output structure
sREC.inSampleRate   = inSampleRate;
sREC.db2LFP         = db2LFP;
sREC.bl1Pres        = bl1Pres;
sREC.bl1Run         = bl1Run;