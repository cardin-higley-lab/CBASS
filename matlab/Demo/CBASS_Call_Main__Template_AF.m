% clear, close all
%% Adds the code to the path ---------- change to the path to the code on your computer -----------------
chCodePath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts';
addpath(genpath(chCodePath)); %Add the root folder to current folder

%% Load the data ----- change to load your data --------------------------------------------------------- 
% The section must output the variable db2LFP (i.e.  a channel x time_sample
% matrix) and its sample rate inSample rate

% Path to the data
% chDataPath      = 'D:\gamma_bouts\';
chDataPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/Data'; % root directory
chExperiment    = 'Example_1';
chExpPath       = fullfile(chDataPath, chExperiment);

% Loads the data
sREC            = CBASS_L0_LoadData(chExpPath);
db2LFP          = sREC.db2LFP;
inSampleRate    = sREC.inSampleRate;
%% Sets the output folder ----- edit or comment to provide externally throug a system call -------------
chOutPath   = fullfile(chDataPath, 'CBASS_Call_Main');
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
[sFREQ_BAND, cTROUGH] = CBASS_Main_DetectEvents(sREC.db2LFP, sREC.inSampleRate, cBAND, cSTATE, sOPTION);

%% Calls embedding and its plotting 

%Address for the output ----- comment to provide externally throug a system call ---------------------
% chOutPath   = fullfile(chPath, 'Figures', 'Figure_4_Embedding');
chTmpPath   = chOutPath;
chOutFile   = chExperiment;

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

% Sets the data format
% chDataFormat    = cFORMAT;
chDataFormat    = 'complex';

% Call python to do the embeding
blZScore        = true;
if blZScore, chZFlag = 'ZScore'; else, chZFlag = 'Raw'; end

% Sets the embeding of choice for
chEmbedMethod   = 'pca';
chMethod = chEmbedMethod;
inN_Component   = 2;
iBnd=2; %1-beta; 2-gamma
chLabel         = sOPTION.cBAND_LABEL{iBnd};

chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat '_' chZFlag];

disp('Embedding...')
[status, commandOut] = CBASS_L2_Embed(cTROUGH{iBnd}, chTmpPath, chOutPath, chExpLabel, ...
    chEmbedMethod, chDataFormat, blZScore, inN_Component)

% Generates the embeding plots for running and not running
disp('Plotting embedding')
in1EmbedLabel   = sREC.bl1Run(cTROUGH{iBnd}.in1Index);
chLabelTag      = ['all_regions_'  sOPTION.cSTATE_LBL{iBnd}];
chFormatImg     = 'png';
chRotate3D      = 'False';
chAddLegend     = 'False';
inFontSize      = '12';
chDiscrete      = 'True';
[status, commandOut] = CBASS_L2_PlotEmbed(cTROUGH{iBnd}, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
    chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);

disp('Plotting embedding')
in1EmbedLabel   = sFREQ_BAND(iBnd).db1Score';
chLabelTag      = ['probTrough_'  sOPTION.cSTATE_LBL{iBnd}];
chFormatImg     = 'png';
chRotate3D      = 'False';
chAddLegend     = 'False';
inFontSize      = '12';
chDiscrete      = 'False';
[status, commandOut] = CBASS_L2_PlotEmbed(cTROUGH{iBnd}, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
    chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);

disp('Plotting embedding')
in1EmbedLabel   = sFREQ_BAND(iBnd).db1Score';
chLabelTag      = ['probTrough_legend_'  sOPTION.cSTATE_LBL{iBnd}];
chFormatImg     = 'eps';
chRotate3D      = 'False';
chAddLegend     = 'True';
inFontSize      = '12';
chDiscrete      = 'False';
[status, commandOut] = CBASS_L2_PlotEmbed(cTROUGH{iBnd}, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
    chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);

disp('Plotting embedding')
in1EmbedLabel   = sFREQ_BAND(iBnd).bl1Partition';
chLabelTag      = ['partition_'  sOPTION.cSTATE_LBL{iBnd}];
chFormatImg     = 'png';
chRotate3D      = 'False';
chAddLegend     = 'False';
inFontSize      = '12';
chDiscrete      = 'False';
[status, commandOut] = CBASS_L2_PlotEmbed(cTROUGH{iBnd}, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
    chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        
%% Saves the ouput structure
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end
save(fullfile(chOutPath, chOutFile), 'sFREQ_BAND', '-v7.3');