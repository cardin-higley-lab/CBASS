%% Script to plot Figure 1 of the paper
% This script exemplifies how to process and plot the Figure 1 (methods).

% clear, close all
%% Loading data ----- change section as needed ----------------------------------------------------------
cEXP = {'Example_1'};

iExp = 1;
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts'; % root directory
% chDir       = fullfile(chPath,'Data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
addpath(genpath(chPath)); %Add the root folder to current folder

% Loads the data
sREC        = CBASS_L0_LoadData(chDir);
sREC        = CBASS_L1_AddPhaseRandomizedSignal(sREC);

%% Address for the output ----- comment to provide externally throug a system call ---------------------
chOutPath   = fullfile(chPath, 'Figures', 'Figure_4_Embedding');
chTmpPath   = chOutPath;
chOutFile   = chExperiment;

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

%% Sets data formatting and embeding global parameters
% Sets the data format
% chDataFormat    = cFORMAT;
chDataFormat    = 'complex';

% Call python to do the embeding
blZScore        = true;
if blZScore, chZFlag = 'ZScore'; else, chZFlag = 'Raw'; end

% Sets the embeding of choice for
chEmbedMethod   = 'umap';
chMethod = chEmbedMethod;
inN_Component   = 3;


%% Set up the loop ----- comment to provide externally throug a system call ---------------------
cBAND_LABEL = {'Beta', 'Gamma'};
cSTATE_LBL  = {'Stim', 'Running'};

% Gets the number of band
inNBnd      = length(cBAND_LABEL);

% Sets reference channels
inRefChan       = 5;

% Sets methods for the construction of the adjacency matrix
cMETHOD     = {'FixedThreshold', 'PerNodeThreshold'};
in1Method   = [1];% 2];

% Sets clustering options
blUseRate   = false;
inNClu      = 20;
%% Calculates filters sets for beta and gamma
for iBnd = 1:inNBnd
    % Sets label
    chLabel         = cBAND_LABEL{iBnd};
    if strcmp(cBAND_LABEL{iBnd},'Beta')
        db1Band       = [15 30];
    elseif strcmp(cBAND_LABEL{iBnd},'Gamma')
        db1Band       = [30 80];
    end
    
    % Sets state
    if strcmp(cSTATE_LBL{iBnd}, 'Stim')
        bl1State      = sREC.bl1Pres;
    elseif strcmp(cSTATE_LBL{iBnd}, 'Running')
        bl1State      = sREC.bl1Run;
    end
    
    % START PROCESSING ------------------------------------------------
    % Computes the trough
    sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
    sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, db1Band, inRefChan, chLabel, chDataFormat);
    
    for iMth = in1Method
        
        %         % Identifies enriched regions
        %         sFILTER     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, cSTATE{iBnd}, blZScore, iMth, inNClu); %Original
        
        [sFILTER, sCLU, in1CluKM, in1Sel] = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, bl1State, blZScore, iMth, inNClu);
        if isempty(sFILTER), sprintf('No Filter for %s %s Method %s', chExperiment, cBAND_LABEL{iBnd}, iMth); continue; end
        inNFlt      = length(sFILTER);
        cFLT_LBL    = cellfun(@(x) sprintf('Filter %d', x), num2cell(1:inNFlt), 'UniformOutput', false);
        
        chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat '_' chZFlag];
        
        disp('Embedding...')
        [status, commandOut] = CBASS_L2_Embed(sTROUGH, chTmpPath, chOutPath, chExpLabel, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component)
        
        %         % Generates the embeding plots for running and not running
        %         disp('Plotting embedding')
        %         in1EmbedLabel   = sREC.bl1Run(sTROUGH.in1Index);
        %         chLabelTag      = 'all_regions';
        %         [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
        %             chMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel);
        
        % Generates the embeding plots for running and not running
        disp('Plotting embedding')
        in1EmbedLabel   = sREC.bl1Run(sTROUGH.in1Index);
        chLabelTag      = ['all_regions_'  cSTATE_LBL{iBnd}];
        chFormatImg     = 'png';
        chRotate3D      = 'False';
        chAddLegend     = 'False';
        inFontSize      = '12';
        chDiscrete      = 'True';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        disp('Plotting embedding')
        ratio_all_regions = NaN(size(sREC.bl1Run(sTROUGH.in1Index)));
        for cluster_num=1:size(sCLU,2)
            ratio_all_regions(sCLU(cluster_num).bl1Member) = sCLU(cluster_num).dbRate;
        end
        in1EmbedLabel   = ratio_all_regions;
        chLabelTag      = ['enriched_'  cSTATE_LBL{iBnd}];
        chFormatImg     = 'eps';
        chRotate3D      = 'False';
        chAddLegend     = 'True';
        inFontSize      = '12';
        chDiscrete      = 'False';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        in1EmbedLabel   = ratio_all_regions;
        chLabelTag      = ['enriched_'  cSTATE_LBL{iBnd}];
        chFormatImg     = 'png';
        chRotate3D      = 'False';
        chAddLegend     = 'False';
        inFontSize      = '12';
        chDiscrete      = 'True';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        in1EmbedLabel   = in1CluKM';
        chLabelTag      = ['20Clusters' cSTATE_LBL{iBnd}];
        chFormatImg     = 'png';
        chRotate3D      = 'False';
        chAddLegend     = 'False';
        inFontSize      = '12';
        chDiscrete      = 'True';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        chFormatImg     = 'eps';
        chRotate3D      = 'False';
        chAddLegend     = 'True';
        inFontSize      = '12';
        chDiscrete      = 'True';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        in1EmbedLabel   = sFILTER.bl1Member';
        chLabelTag      = ['SigClusters' cSTATE_LBL{iBnd}];
        chFormatImg     = 'eps';
        chRotate3D      = 'False';
        chAddLegend     = 'True';
        inFontSize      = '12';
        chDiscrete      = 'True';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
        
        chFormatImg     = 'png';
        chRotate3D      = 'False';
        chAddLegend     = 'False';
        inFontSize      = '12';
        chDiscrete      = 'True';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete);
    end
    close all
end