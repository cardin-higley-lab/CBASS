%% Loads data
clear, close all
chExperiment = 'Example_1';
chPath      = '/gpfs/ysm/home/ahf38/project/gamma_bouts/Python_Sanbox/';
chDir       = fullfile(chPath,'data', chExperiment);
% chPath      = 'D:\gamma_bouts\';
% chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Tests', 'Test_PolarVsComplex_RawVsZScored');
chTmpPath   = fullfile(chPath,'Tests', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

%% Sets embeding global parameters
chMethod        = 'phate';
inN_Component   = 3;
%% Loads the data
sREC        = CBASS_L0_LoadData(chDir);
sREC        = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Set up the loop
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND       = {[15 30], [30 80]};
cSTATE      = {sREC.bl1Pres sREC.bl1Run};
cSTATE_LBL  = {'Stim', 'Running'};
cFORMAT     = {'complex', 'polar', 'complex', 'polar'};
bl1ZScore   = [false false true true];

%% Loops through conditons
for iBnd = 1:length(cBAND)
    for iFmt = 1:length(cFORMAT)
        chLabel         = cBAND_LABEL{iBnd};
        dbBand          = cBAND{iBnd};
        inRefChan       = 5;
        chDataFormat    = cFORMAT{iFmt};
        
        % Call python to do the embeding
        blZScore        = bl1ZScore(iFmt);
        if blZScore, chZFlag = 'ZScore'; else, chZFlag = 'Raw'; end
        
        % Sets the gloabal label for the experiment
        chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat '_' chZFlag];
        
        % START PROCESSING ------------------------------------------------
        % Computes the trough
        sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
%         sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
        
        [status, commandOut] = CBASS_L2_Embed(sTROUGH, chTmpPath, chOutPath, chExpLabel, ...
            chMethod, chDataFormat, blZScore, inN_Component)
        
        % Identifies enriched regions
        [sREGION, sCLU, sFILTERS, cBRANCH] = CBASS_L2_GetEnrichedRegion(sTROUGH, sREC.db2LFP, sREC.inSampleRate, cSTATE{iBnd}, blZScore);
        if isempty(sREGION); fprintf('No enriched region for %s cond %d\r', cBAND_LABEL{iBnd}, iFmt); continue; end
        
        % Generates the embeding plots
        in1EmbedLabel   = sREGION(1).bl1Member';
        chLabelTag      = 'Region1';
        [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
            chMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel);
        
        % Get the pulses
        for iREGIONS = 1:length(sREGION)
%             iREGIONS
            sPULSE{iREGIONS}      = CBASS_L3_GetPulse(sREC, sTROUGH, sREGION(iREGIONS).bl1Member);
        end
        
        % PLOT AND SAVE FIGURES -------------------------------------------
%         % 1 -- Plots the histogram for the pulses
%         hFIG = figure;
%         chFigName = ['Peak_Histogram_' chExpLabel];
%         CBASS_Plot_PeakHistogram(sPULSE);
%         CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
%         % 2 -- Sets the event of interst and plots event triggered average
%         cEVENT          = {sPULSE.bl1Pulse sTROUGH.in1Index sTROUGH.in1Index(sREGION(1).bl1Member)};
%         cEVENT_LABEL    = {'Pulse', 'Everything', 'SelTrough'};
%         inNEvt          = length(cEVENT);
%         
%         hFIG = figure('Position', [50 50 1250 600]);
%         chFigName = ['Triggered_Averages_' chExpLabel];
%         for iEvt = 1:inNEvt
%             subplot(2, inNEvt, iEvt)
%             CBASS_Plot_LFP_EventTriggeredAverage(sREC.db2LFP, cEVENT{iEvt}, [-.075 .075], sREC.inSampleRate);
%             title(sprintf('%s Triggered Average: LFP', cEVENT_LABEL{iEvt}))
%             subplot(2, inNEvt, inNEvt + iEvt)
%             CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, cEVENT{iEvt}, [-.075 .075], sREC.inSampleRate);
%             title('CSD')
%         end
%         CBASS_SaveFig(chOutPath, hFIG, {chFigName});
%         
        % 2.1 -- Sets the event of interst and plots event triggered
        % average for each pulse
%         cEVENT          = sPULSE.bl1Pulse;
        cEVENT_LABEL    = {'Pulse1', 'Pulse2', 'Pulse3', 'Pulse4', 'Pulse5', 'Everything'};
        inNEvt          = length(sPULSE)+1;
        
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = ['Triggered_Averages_' chExpLabel];
        for iEvt = 1:inNEvt
            subplot(2, inNEvt, iEvt)
            try
                CBASS_Plot_LFP_EventTriggeredAverage(sREC.db2LFP, sPULSE{iEvt}.bl1Pulse, [-.075 .075], sREC.inSampleRate);
            catch
                CBASS_Plot_LFP_EventTriggeredAverage(sREC.db2LFP, sTROUGH.in1Index, [-.075 .075], sREC.inSampleRate);
            end
            title(sprintf('%s Triggered Average: LFP', cEVENT_LABEL{iEvt}))
            subplot(2, inNEvt, inNEvt + iEvt)
            try
                db2CSD_ETA_Inter{iEvt} = CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sPULSE{iEvt}.bl1Pulse, [-.075 .075], sREC.inSampleRate);
            catch
                db2CSD_ETA_Inter{iEvt} = CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, sTROUGH.in1Index, [-.075 .075], sREC.inSampleRate);
            end
            title('CSD')
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
        % 2.2 -- Compute the correlation between the pulses
%         data = nan(length(sFILTERS),numel(sFILTERS{i}));
%         inNEvt = length(sFILTERS);
%         ConfusionMatrix = nan(inNEvt,inNEvt);
%         for i=1:inNEvt
%             for j=1:inNEvt
% %                 if i==j
% %                     ConfusionMatrix(i,j)= 0;
% %                 else
%                 aux = corrcoef(sFILTERS{i}, sFILTERS{j},'rows','complete');
%                 ConfusionMatrix(i,j)= aux(1,2);
% %                 end
% %                 ConfusionMatrix(j,i)= aux(1,2);
%                 data(i,:) = sFILTERS{i}(:);
%             end
%         end
%         figure, heatmap(ConfusionMatrix)
%         data = ConfusionMatrix;
%         %% Affinity matrix
%         sigma=0.05
%         for i=1:size(data,1)    
%             for j=1:size(data,1)
%                 if i==j
%                     affinity(i,j) = 0;
%                 else
%                     dist = sqrt((data(i,1) - data(j,1))^2 + (data(i,2) - data(j,2))^2); 
%                     affinity(i,j) = exp(-dist/(2*sigma^2));
%                 end
%             end
%         end
        % 3 -- Plot the state enrichment 
        cEVENT          = {sPULSE.bl1Pulse sTROUGH.in1Index(sREGION(1).bl1Member)};
        cEVENT_LABEL    = {'Pulse', 'SelTrough'};
        inNEvt          = length(cEVENT);
        
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = ['PulseVs' cSTATE_LBL{iBnd} '_' chExpLabel];
        for iEvt = 1:inNEvt
            subplot(1, 2, iEvt)
            chSigString = CBASS_Plot_Pulse_vs_State(cEVENT{iEvt}, cSTATE{iBnd}, [-5 10], sREC.inSampleRate);
            title(sprintf('%s : %s', cEVENT_LABEL{iEvt}, chSigString))
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
        % 4 -- Plots LFP example
        in1Anchor       = find(cSTATE{iBnd}(1:end - 1) & cSTATE{iBnd}(2:end));
        cEVENT          = {sPULSE.bl1Pulse  sTROUGH.in1Index(sREGION(1).bl1Member)};
        cEVENT_LABEL    = {'Pulse', 'SelTrough'};
        
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = ['LFP_Examples_' chExpLabel];
        for iEx = 1:4
            subplot(2, 2, iEx)
            CBASS_Plot_LFP_EventExample(sREC.db2LFP, in1Anchor(iEx), [-2 3], sREC.inSampleRate, cEVENT, cEVENT_LABEL)
            title(sprintf('Data around %s onset: %d', cSTATE_LBL{iBnd}, iEx));
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
        close all
    end
end
