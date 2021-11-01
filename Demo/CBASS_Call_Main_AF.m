% clear, close all
%% Loading data ----- change section as needed ---------------------------------------------------------- 
cEXP = {'Example_1'}; %, 'Example_2'}; % List of experiments 

iExp = 1;
chExperiment = cEXP{iExp};
chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts'; % root directory
chDir       = fullfile(chPath,'Data', chExperiment); 
% chPath      = 'D:\gamma_bouts\';
% chDir       = fullfile(chPath, chExperiment);

addpath(genpath(chPath)); %Add the root folder to current folder

% Loads the data
sREC            = CBASS_L0_LoadData(chDir);

%% Address for the output ----- comment to provide externally throug a system call ---------------------
chOutPath   = fullfile(chDir, 'CBASS_Call_Main');
chOutFile   = chExperiment;

%% Non optional parameters ---- comment to provide externally throug a system call ---------------------
db2LFP          = sREC.db2LFP;
inSampleRate    = sREC.inSampleRate;
% cBAND           = {[15 30], [30 80]};
% cSTATE          = {sREC.bl1Pres sREC.bl1Run};

%% Optional parameters ----- Uncomment to provide externally through a system call using carrier name -- 
sOPTION.cBAND_LABEL     = cBAND_LABEL;  % character array or cell array of labels  for bands in cBAND (i.e. {'Beta', 'Gamma'})
sOPTION.cSTATE_LBL      = cSTATE_LBL  % character array or cell array of labels  for states in cSTATE (i.e. {'Stim', 'Running'})
% sOPTION.blVerbose       = true;    % Sets whether to print progress on the screen. Default is true 

% L1 options for formatting hilbert troughs (Function CBASS_L1_GetTrough)
sOPTION.chDataFormat    = cFORMAT; % Format of hilbert transform coordinates ('polar' or 'complex' Default is 'complex')
% sOPTION.inRefChan       = 5;    % Reference channel (Default is the last of row of db2LFP.

% L2 options for the computation of filters (Function CBASS_L2_GetFilters)
% sOPTION.blZScore        = bl1ZScore;     % z-score for k-means partitioning. Default is true.
% sOPTION.inMethod        = 1;     % method for spectral clusering of template. Default is 1.
% sOPTION.inNClu          = 20;       % number of cluster of k-means partitioning. Default is 20.
% sOPTION.dbSigThrs       = 10.^-4;    % threshold for enrichment significance. Default is 10.^-4.
% sOPTION.inNMaxIter      = 10000;   % maximum iteration of k-means partioning. Default is 10000;

% L3 options for template matching
% sOPTION.blCntr          = true;       % Sets whether to center signal and template about their mean. Default is true
% sOPTION.blNorm          = true;       % Sets whether to normalize signal and template by their S.D. Default is true.

if exist('sOPTION','var')
    chFieldNames = struct2cell(sOPTION);
    for idxField=1:length(chFieldNames)
        chOutFile = [chOutFile '_' chFieldNames{idxField}];
    end
end

% Setting the cSTATE variable according to the external definition
if isfield(sOPTION,'cSTATE_LBL')
    if strcmp(sOPTION.cSTATE_LBL,'Stim')
        cSTATE  = {sREC.bl1Pres};
    elseif strcmp(sOPTION.cSTATE_LBL,'Running')
        cSTATE = {sREC.bl1Run};
    end
end
disp(['File name: ' chOutFile])
%% Calls main function ----- Add sOPTION if any line above was uncommented  --
sFREQ_BAND = CBASS_Main_DetectEnrichedMotif(sREC.db2LFP, sREC.inSampleRate, cBAND, cSTATE);

%% Saves the ouput structure
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end
save(fullfile(chOutPath, chOutFile), 'sFREQ_BAND', '-v7.3');