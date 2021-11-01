%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
for iExp = 1:length(cEXP)
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Tests', 'Test_Clustering_Method_2');
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


bl1UseRate  = [false(1, 12) true(1, 12)];
in1Method   = [1 1 1 1 2 2 2 2 3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3];
in1NClu     = repmat([10 20 40 80], 1, 6);
inNFmt      = length(bl1UseRate);

[cFORMAT, cFMT_LBL] = deal(cell(1, inNFmt));
for iFmt = 1:inNFmt
    if bl1UseRate(iFmt), chS1 = 'Rate'; chS2 = 'R'; else, chS1 = 'NoRate'; chS2 = 'NR'; end
    cFORMAT{iFmt}   = sprintf('v%d_%dClu_%s', in1Method(iFmt), in1NClu(iFmt), chS1);
    cFMT_LBL{iFmt}  = sprintf('%dCv%d%s', in1NClu(iFmt), in1Method(iFmt), chS2);
end
%% Loops through conditons
for iBnd = 1:length(cBAND)
% for iBnd = 2
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
    cEVENT = {};
    
    % Loops over formats
    for iFmt = 1:inNFmt
        
        % Sets option
        blUseRate   = bl1UseRate(iFmt);
        inMethod    = in1Method(iFmt);
        inNClu      = in1NClu(iFmt);
   
        % Identifies enriched regions
        [sREGION, sCLU] = CBASS_L2_GetEnrichedRegion(sTROUGH, cSTATE{iBnd}, blZScore, blUseRate, inMethod, inNClu);
        if isempty(sREGION); fprintf('No enriched region for %s cond %d\r', cBAND_LABEL{iBnd}, iFmt); continue; end
        
        % Generates the embeding plots
        %         in1EmbedLabel   = sREGION(1).bl1Member';
        %         chLabelTag      = 'Region1';
        %         [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExpLabel, chLabelTag, ...
        %             chEmbedMethod, chDataFormat, blZScore, inN_Component, in1EmbedLabel)
        
        % Get the pulses for each enriched region
        sPULSE = CBASS_L3_GetPulse(sREC, sTROUGH, sREGION(1).bl1Member);
        
        % Aggregates the pulse events with the most enriched region
        cEVENT = cat(2, cEVENT, {sPULSE.bl1Pulse});
    end
    
    inNFig = ceil(inNFmt / 8);
    for iFig = 1:inNFig
        hFIG = figure('Position', [50 50 1250 600]);
        chFigName = ['CSD_' chExpLblShrt '_' num2str(iFig)];
        for iPlt = 1:8
            iFmt = iPlt + (8 * (iFig - 1));
            subplot(2, 4, iPlt)
            CBASS_Plot_CSD_EventTriggeredAverage(sREC.db2LFP, cEVENT{iFmt}, [-.075 .075], sREC.inSampleRate);
            title(sprintf('%s : CSD', cFMT_LBL{iFmt}))
            if iFmt >= inNFmt; break; end
        end
        CBASS_SaveFig(chOutPath, hFIG, {chFigName});
    end
    
    % Plots the first figure
    clear hPLT
    hFIG = figure('Position', [50 50 1250 600]); hold on
    chFigName = ['Enrichment_' cSTATE_LBL{iBnd} '_' chExpLblShrt];
    [db1FreqOFF, db1FreqON, db1Sig] = deal(nan(1, inNFmt));
    for iFmt = 1:inNFmt
        % Frequency
        db1FreqON(iFmt)    = sREC.inSampleRate * sum(cEVENT{iFmt} & cSTATE{iBnd})./sum(cSTATE{iBnd});
        db1FreqOFF(iFmt)   = sREC.inSampleRate * sum(cEVENT{iFmt} & ~cSTATE{iBnd})./sum(~cSTATE{iBnd});
        
        % Significance of the increase
        dbPVal = 2 * (1 - normcdf(abs(db1FreqOFF(iFmt) - db1FreqON(iFmt))./...
            sqrt((db1FreqOFF(iFmt) ./ sum(cEVENT{iFmt} & cSTATE{iBnd})) +  ...
            (db1FreqON(iFmt) ./ sum(cEVENT{iFmt} & ~cSTATE{iBnd})) )));
        if dbPVal < .05; db1Sig(iFmt) = 1; end
    end
    hPLT(1) = plot(1:inNFmt, db1FreqOFF, 'bo');
    hPLT(2) = plot(1:inNFmt, db1FreqON, 'rx');
    dbLevel = max([db1FreqON db1FreqOFF]) + .1 * range([db1FreqON db1FreqOFF]);
    plot(1:inNFmt, db1Sig * dbLevel, 'k*')
    ylabel('Frequency (Hz)'), xlabel('Condition'); title('Pulse Frequency');
    xlim([0 inNFmt + 1]), set(gca, 'XTick', 1:inNFmt, 'XTickLabel', cFMT_LBL);
    xtickangle(90);
    legend(hPLT, {[cSTATE_LBL{iBnd} ' OFF'], [cSTATE_LBL{iBnd} ' ON']})
    CBASS_SaveFig(chOutPath, hFIG, {chFigName});
    
    % Plots the second figure
    hFIG = figure('Position', [50 50 1250 600]);
    chFigName = ['Overlap_' cSTATE_LBL{iBnd} '_' chExpLblShrt];
    cEVENT_IDX = cellfun(@find, cEVENT, 'UniformOutput', false);
    db2D_Asym = EventOverlap(cEVENT_IDX, .005, sREC.inSampleRate);
    imagesc(1:inNFmt, 1:inNFmt, db2D_Asym)
    ylabel('Cond'); xlabel('Cond'); title('EventOverlap');
    set(gca, 'XTick', 1:inNFmt, 'XTickLabel', cFMT_LBL, 'YTick', 1:inNFmt, 'YTickLabel', cFMT_LBL);
    xtickangle(90);
    colorbar
    CBASS_SaveFig(chOutPath, hFIG, {chFigName});
    
    close all
end
end