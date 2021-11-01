%% Loads data
clear, close all
cEXP = {'Example_1'};%, 'Example_2'};
for iExp = 1:length(cEXP)
% iExp = 1;
chExperiment = cEXP{iExp};
chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts'; % root directory
chDir       = fullfile(chPath,'Data', chExperiment);
% chPath      = 'D:\gamma_bouts\';
% chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Pipeline', 'Test_Filters_Conv');
chTmpPath   = fullfile(chPath,'Pipeline', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder
addpath(fullfile(chPath, 'Pipeline', 'Utilities'))

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

%% Sets data formatting and embeding global parameters
% Sets the data format
chDataFormat    = 'complex';

% Call python to do the embeding
blZScore        = true;
if blZScore, chZFlag = 'ZScore'; else, chZFlag = 'Raw'; end

% Sets the embeding of choice for
chEmbedMethod   = 'umap';
chMethod = chEmbedMethod;
inN_Component   = 3;
%% Loads the data
sREC        = CBASS_L0_LoadData(chDir);
sREC        = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Set up the loop
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND       = {[15 30], [30 80]};
cSTATE      = {sREC.bl1Pres sREC.bl1Run};
cSTATE_LBL  = {'Stim', 'Running'};
inNBnd      = length(cBAND_LABEL);

% Sets methods for the construction of the adjacency matrix
cMETHOD     = {'FixedThreshold', 'PerNodeThreshold'};
in1Method   = [1];% 2];

% Sets clustering options
blUseRate   = false;
inNClu      = 20;
%% Calculates filters sets for beta and gamma
for iBnd = 1:inNBnd
    chLabel         = cBAND_LABEL{iBnd};
    dbBand          = cBAND{iBnd};
    inRefChan       = 5;
    
    % START PROCESSING ------------------------------------------------
    % Computes the trough
    sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
    sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
    
    for iMth = in1Method

        % Identifies enriched regions
%         sFILTER     = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, cSTATE{iBnd}, blZScore, iMth, inNClu); %Original 


        
        [sFILTER, sCLU, in1CluKM, in1Sel] = CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, cSTATE{iBnd}, blZScore, iMth, inNClu);
        if isempty(sFILTER), sprintf('No Filter for %s %s Method %s', chExperiment, cBAND_LABEL{iBnd}, iMth); continue; end
        inNFlt      = length(sFILTER);
        cFLT_LBL    = cellfun(@(x) sprintf('Filter %d', x), num2cell(1:inNFlt), 'UniformOutput', false);
        
        chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat '_' chZFlag];
         
        disp('Embedding...')
        [status, commandOut] = CBASS_L2_Embed(sTROUGH, chTmpPath, chOutPath, chExpLabel, ...
        chEmbedMethod, chDataFormat, blZScore, inN_Component)
    
        % Generates the embeding plots for running and not running
        disp('Plotting embedding')
        in1EmbedLabel   = sREC.bl1Run(sTROUGH.in1Index);
        chLabelTag      = 'all_regions';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel);
        
        % Generates the embeding plots
        disp('Plotting embedding')
        in1EmbedLabel   = sREGION(1).bl1Member';
        chLabelTag      = 'Region1';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel);
        
        % Plots the filters
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = [cBAND_LABEL{iBnd} '_' num2str(iMth) '_' chExperiment '_Filters'];
        [inRow, inCol] = FindPlotNum(length(sFILTER));
        for iFlt = 1:inNFlt
            subplot(inRow, inCol, iFlt)
            imagesc(sFILTER(iFlt).db2Filter);
            xlabel(strrep(cFLT_LBL{iFlt}, '_', '\_'))
            if iFlt == round(inCol/2); title(strrep(chFigName, '_', '\_')); end
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
        % Gets the pulse
        clear sPULSE
        for iFlt = 1:inNFlt
            % Gets the pulse
            sPULSE(iFlt) = CBASS_L3_GetPulse_2(sREC, sFILTER(iFlt).db2Filter);
        end
        
        if iFlt > 1
            % Finds the overlap between the pulse
            cEVENT_IDX = cellfun(@find, {sPULSE.bl1Pulse}, 'UniformOutput', false);
            [db2D_Asym] = EventOverlap(cEVENT_IDX, .005, sREC.inSampleRate);
            
            hFIG = figure('Position', [50 50 1250 600]);
            chFigName = [cBAND_LABEL{iBnd} '_' num2str(iMth) '_' chExperiment '_PulseOverlap'];
            imagesc(1:inNFlt, 1:inNFlt, db2D_Asym)
            ylabel('Filter'); xlabel('Filter'); title('EventOverlap');
            xtickangle(90);
            colorbar
            CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        end
        
        % Plots the Pulse vs state
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = [cBAND_LABEL{iBnd} '_' num2str(iMth) '_' chExperiment '_PulseVsState'];
        for iFlt = 1:inNFlt
            subplot(inRow, inCol, iFlt)
            chSigString = CBASS_Plot_Pulse_vs_State(sPULSE(iFlt).bl1Pulse, cSTATE{iBnd}, [-5 10], sREC.inSampleRate);
            title(sprintf('%s : %s', cFLT_LBL{iFlt}, chSigString))
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
    end
    close all
end
end