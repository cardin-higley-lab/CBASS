%% Loads data
clear, close all
chExperiment = 'Example_1';
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Tests', 'Test_Clustering_Method');
chTmpPath   = fullfile(chPath,'Tests', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder

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
inN_Component   = 3;
%% Loads the data
sREC        = CBASS_L0_LoadData(chDir);
sREC        = CBASS_L1_AddPhaseRandomizedSignal(sREC);
%% Set up the loop
cBAND_LABEL = {'Beta', 'Gamma'};
cBAND       = {[15 30], [30 80]};
cSTATE      = {sREC.bl1Pres sREC.bl1Run};
cSTATE_LBL  = {'Stim', 'Running'};

cFORMAT     = {'v1_NoRate', 'v2_NoRate', 'v3_NoRate', 'v1_Rate', 'v2_Rate', 'v3_Rate'};
cFMT_LBL    = {'v1NR', 'v2NR', 'v3R', 'v1R', 'v2R', 'v3R'};
bl1UseRate  = [false false false true true true];
in1Method   = [1 2 3 1 2 3];
inNFmt      = length(cFORMAT);
%% Loops through conditons
for iBnd = 1:length(cBAND)
    chLabel         = cBAND_LABEL{iBnd};
    dbBand          = cBAND{iBnd};
    inRefChan       = 5;
    
    % Sets the gloabal label for the experiment
    chExpLabel      = [chExperiment '_' chLabel '_' chDataFormat '_' chZFlag];
    chExpLblShrt    = [chExperiment '_' chLabel];
    
    % START PROCESSING ------------------------------------------------
    % Computes the trough
    sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, sREC.inSampleRate, dbBand, inRefChan, chLabel, chDataFormat);
    
    % Performs the embeding if the emmbeding file does not exist
    %         [status, commandOut] = CBASS_L2_Embed(sTROUGH, chTmpPath, chOutPath, chExpLabel, ...
    %             chMethod, chDataFormat, blZScore, inN_Component)
    
    % Initizalize aggregate variable
    cEVENT_1 = {};
    
    % Loops over formats
    for iFmt = 1:inNFmt
        
        % Sets option
        blUseRate   = bl1UseRate(iFmt);
        inMethod    = in1Method(iFmt);
        
        % Identifies enriched regions
        [sREGION, sCLU] = CBASS_L2_GetEnrichedRegion(sTROUGH, cSTATE{iBnd}, blZScore, blUseRate, inMethod);
        if isempty(sREGION); fprintf('No enriched region for %s cond %d\r', cBAND_LABEL{iBnd}, iFmt); continue; end
        
        % Generates the embeding plots
        %         in1EmbedLabel   = sREGION(1).bl1Member';
        %         chLabelTag      = 'Region1';
        %         [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
        %             chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel)
        
        % Get the pulses for each enriched region
        inNReg = length(sREGION);
        clear sPULSE
        for iReg = 1:inNReg
            sPULSE(iReg)    = CBASS_L3_GetPulse(sREC, sTROUGH, sREGION(iReg).bl1Member);
        end
        
        % PLOT AND SAVE FIGURES -------------------------------------------
        % 1 -- Plots the histogram for the pulses
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = ['Peak_Histogram_' chExpLblShrt '_' cFORMAT{iFmt}];
        [inNRow, inNCol] = FindPlotNum(inNReg);
        for iReg = 1:inNReg
            subplot(inNRow, inNCol, iReg);
            CBASS_Plot_PeakHistogram(sPULSE(iReg), iReg == 1);
            title(sprintf('Region %d', iReg))
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
        % 2 -- Sets the event of interst and plots event triggered average
        cEVENT          = {sPULSE.bl1Pulse};
        cEVENT_LABEL    = cellfun(@(x) sprintf('Region %d', x), num2cell(1:inNReg), 'UniformOutput', false);
        inNFig = ceil(inNReg / 4);
        for iFig = 1:inNFig
            hFIG = figure('Position', [50 50 1250 600]);
            chFigName = ['Triggered_Averages_' chExpLblShrt '_' cFORMAT{iFmt} '_' num2str(iFig)];
            for iPlt = 1:4
                iReg = iPlt + (4 * (iFig - 1));
                subplot(2, 4, iPlt)
                CBASS_Plot_LFP_EventTriggeredAverage(sREC.db2LFP, cEVENT{iReg}, [-.075 .075], sREC.inSampleRate);
                title(sprintf('%s Triggered Average: LFP', cEVENT_LABEL{iReg}))
                subplot(2, 4, iPlt + 4)
                CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, cEVENT{iReg}, [-.075 .075], sREC.inSampleRate);
                title('CSD')
                if iReg >= inNReg; break; end
            end
            CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        end
        
        clear hPLT
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = ['Overlap_And_Enrichment_' cSTATE_LBL{iBnd} '_' chExpLblShrt '_' cFORMAT{iFmt}];
        [db1FreqOFF, db1FreqON, db1Sig] = deal(nan(1, inNReg));
        for iReg = 1:inNReg
            % Frequency
            db1FreqON(iReg)    = sREC.inSampleRate * sum(cEVENT{iReg} & cSTATE{iBnd})./sum(cSTATE{iBnd});
            db1FreqOFF(iReg)   = sREC.inSampleRate * sum(cEVENT{iReg} & ~cSTATE{iBnd})./sum(~cSTATE{iBnd});
            
            % Significance of the increase
            dbPVal = 2 * (1 - normcdf(abs(db1FreqOFF(iReg) - db1FreqON(iReg))./...
                sqrt((db1FreqOFF(iReg) ./ sum(cEVENT{iReg} & cSTATE{iBnd})) +  ...
                (db1FreqON(iReg) ./ sum(cEVENT{iReg} & ~cSTATE{iBnd})) )));
            if dbPVal < .05; db1Sig(iReg) = 1; end
        end
        subplot(1, 2, 1); hold on
        hPLT(1) = plot(1:inNReg, db1FreqOFF, 'bo');
        hPLT(2) = plot(1:inNReg, db1FreqON, 'rx');
        dbLevel = max([db1FreqON db1FreqOFF]) + .1 * range([db1FreqON db1FreqOFF]);
        plot(1:inNReg, db1Sig * dbLevel, 'k*')
        ylabel('Frequency (Hz)'), xlabel('Region'); title('Pulse Frequency');
        xlim([0 inNReg + 1])
        legend(hPLT, {[cSTATE_LBL{iBnd} ' OFF'], [cSTATE_LBL{iBnd} ' ON']})
        %CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
        subplot(1, 2, 2)
        cEVENT_IDX = cellfun(@find, cEVENT, 'UniformOutput', false);
        db2D_Asym = EventOverlap(cEVENT_IDX, .005, sREC.inSampleRate);
        imagesc(1:inNReg, 1:inNReg, db2D_Asym)
        ylabel('Region'); xlabel('Region'); title('EventOverlap');
        colorbar
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
        
        % Close existing figures
        close all
        
        % Aggregates the pulse events with the most enriched region
        cEVENT_1 = cat(2, cEVENT_1, cEVENT(1));
    end
    
    clear hPLT
    hFIG = figure('Position', [50 50 1250 600]);
    chFigName = ['Overlap_And_Enrichment_' cSTATE_LBL{iBnd} '_' chExpLblShrt '_Comp'];
    [db1FreqOFF, db1FreqON, db1Sig] = deal(nan(1, inNFmt));
    for iFmt = 1:inNFmt
        % Frequency
        db1FreqON(iFmt)    = sREC.inSampleRate * sum(cEVENT_1{iFmt} & cSTATE{iBnd})./sum(cSTATE{iBnd});
        db1FreqOFF(iFmt)   = sREC.inSampleRate * sum(cEVENT_1{iFmt} & ~cSTATE{iBnd})./sum(~cSTATE{iBnd});
        
        % Significance of the increase
        dbPVal = 2 * (1 - normcdf(abs(db1FreqOFF(iFmt) - db1FreqON(iFmt))./...
            sqrt((db1FreqOFF(iFmt) ./ sum(cEVENT_1{iFmt} & cSTATE{iBnd})) +  ...
            (db1FreqON(iFmt) ./ sum(cEVENT_1{iFmt} & ~cSTATE{iBnd})) )));
        if dbPVal < .05; db1Sig(iFmt) = 1; end
    end
    subplot(1, 2, 1); hold on
    hPLT(1) = plot(1:inNFmt, db1FreqOFF, 'bo');
    hPLT(2) = plot(1:inNFmt, db1FreqON, 'rx');
    dbLevel = max([db1FreqON db1FreqOFF]) + .1 * range([db1FreqON db1FreqOFF]);
    plot(1:inNFmt, db1Sig * dbLevel, 'k*')
    ylabel('Frequency (Hz)'), xlabel('Condition'); title('Pulse Frequency');
    xlim([0 inNFmt + 1]), set(gca, 'XTick', 1:inNFmt, 'XTickLabel', cFORMAT);
    legend(hPLT, {[cSTATE_LBL{iBnd} ' OFF'], [cSTATE_LBL{iBnd} ' ON']})
    %CBASS_SaveFig(chOutPath, hFIG, {chFigName}, 'png');
    subplot(1, 2, 2)
    cEVENT_IDX = cellfun(@find, cEVENT_1, 'UniformOutput', false);
    db2D_Asym = EventOverlap(cEVENT_IDX, .005, sREC.inSampleRate);
    imagesc(1:inNFmt, 1:inNFmt, db2D_Asym)
    ylabel('Cond'); xlabel('Cond'); title('EventOverlap');
    set(gca, 'XTick', 1:inNFmt, 'XTickLabel', cFMT_LBL, 'YTick', 1:inNFmt, 'YTickLabel', cFMT_LBL);
    colorbar
    CBASS_SaveFig(chOutPath, hFIG, {chFigName});
    
    close all
end