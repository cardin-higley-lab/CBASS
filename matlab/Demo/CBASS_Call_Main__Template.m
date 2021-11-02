% clear, close all
%% Adds the code to the path ---------- change to the path to the code on your computer -----------------
chCodePath      = 'D:\gamma_bouts\Pipeline';
addpath(genpath(chCodePath)); %Add the root folder to current folder

%% Load the data ----- change to load your data --------------------------------------------------------- 
% The section must output the variable db2LFP (i.e.  a channel x time_sample
% matrix) and its sample rate inSample rate

% Path to the data
chDataPath      = 'D:\gamma_bouts\';
chExperiment    = 'Example_1';
chExpPath       = fullfile(chDataPath, chExperiment);

% Loads the data
sREC            = CBASS_L0_LoadData(chExpPath);
db2LFP          = sREC.db2LFP;
inSampleRate    = sREC.inSampleRate;
%% Sets the output folder ----- edit or comment to provide externally throug a system call -------------
chOutPath   = fullfile('D:\gamma_bouts\', 'CBASS_Call_Main');
chOutFile   = chExperiment;
%% Non optional parameters ----- edit or comment to provide externally throug a system call -------------
cBAND           = {[15 30], [30 80]};
cSTATE          = {sREC.bl1Pres sREC.bl1Run};
%% Optional parameters ----- Uncomment to provide externally through a system call using carrier name -- 
sOPTION.cBAND_LABEL     = {'Beta', 'Gamma'};  % character array or cell array of labels  for bands in cBAND (i.e. {'Beta', 'Gamma'})
sOPTION.cSTATE_LBL      = {'Stim', 'Running'};  % character array or cell array of labels  for states in cSTATE (i.e. {'Stim', 'Running'})
% sOPTION.blVerbose       = true;    % Sets whether to print progress on the screen. Default is true 

% L1 options for formatting hilbert troughs (Function CBASS_L1_GetTrough)
% sOPTION.chDataFormat    = 'complex'; % Format of hilbert transform coordinates ('polar' or 'complex' Default is 'complex')
% sOPTION.inRefChan       = 5;    % Reference channel (Default is the last of row of db2LFP.

% L2 options for the computation of filters (Function CBASS_L2_GetFilters)
% sOPTION.cBASELINE       = {~sREC.b1lPres ~sREC.bl1Run}; % Baseline state.    
% sOPTION.blZScore        = true;     % z-score for k-means partitioning. Default is true.
% sOPTION.inMethod        = 1;     % method for spectral clusering of template. Default is 1.
% sOPTION.inNClu          = 20;       % number of cluster of k-means partitioning. Default is 20.
% sOPTION.dbSigThrs       = 10.^-4;    % threshold for enrichment significance. Default is 10.^-4.
% sOPTION.inNIter         = 1000;   % maximum iteration of k-means partioning. Default is 10000;
%% Calls main function 
if ~exist('sOPTION', 'var'), sOPTION = struct; end % Creates sOPTION if it does not exist
sFREQ_BAND = CBASS_Main_DetectEvents(sREC.db2LFP, sREC.inSampleRate, cBAND, cSTATE, sOPTION);
%% Saves the ouput structure
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end
save(fullfile(chOutPath, chOutFile), 'sFREQ_BAND', '-v7.3');