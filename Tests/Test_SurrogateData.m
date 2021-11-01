%% Loads data
clear, close all
cEXP = {'Example_1', 'Example_2'};
for iExp = 1:length(cEXP)
chExperiment = cEXP{iExp};
% chPath      = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/';
% chDir       = fullfile(chPath,'data', chExperiment);
chPath      = 'D:\gamma_bouts\';
chDir       = fullfile(chPath, chExperiment);
chOutPath   = fullfile(chPath,'Tests', 'Test_SurrogateData');
chTmpPath   = fullfile(chPath,'Tests', 'tmp');
addpath(genpath(chPath)); %Add the root folder to current folder

%Create output folder if doens't exist
if ~exist(chTmpPath, 'dir'), mkdir(chTmpPath), end
if ~exist(chOutPath, 'dir'), mkdir(chOutPath), end

chDataFormat    = 'complex';

% Call python to do the embeding
blZScore        = true;
if blZScore, chZFlag = 'ZScore'; else, chZFlag = 'Raw'; end

% Sets the embeding of choice for
chEmbedMethod   = 'umap';
inN_Component   = 3;

% Loads the data
sREC        = CBASS_L0_LoadData(chDir);
%% Test the method of surrogate data of different size.
% Set LFP and surrogate signal length
inLenInSec  = 10;
inLenOutSec = 30 * inLenInSec;
inNIn       = (inLenInSec * sREC.inSampleRate);
inNOut      = (inLenOutSec * sREC.inSampleRate);

% Extract a test chunk of LFP
db2LFP_In   = sREC.db2LFP(:, 1:inNIn);
inNChan     = size(db2LFP_In, 1);

% Definese method tag for figures
cMETHOD = {'Interpolated', 'ColoredNoiseFit'};

% Loops through methods
for iMth = 1:length(cMETHOD)
    % Decomposes the LFP into orthogonal principal components
    [db2Coeff, db2LFP_Rnd] = pca(db2LFP_In');
    
    % Randomizes the phase of the resulting signal
    if iMth == 1
        db2LFP_Rnd = CBASS_U_PowerMatchingRandSeries(db2LFP_Rnd, 1, inNOut);
    elseif iMth == 2
        % Alternatively fits a 1/f^alpha model to the spectrum of the referecen
        % data and generate a phase randomized time series
        % Max frequency for the fit;
        dbFreqMax   = 120;
        inNSampFit  = round(120 .* size(db2LFP_In, 2) ./ sREC.inSampleRate);
        [db2LFP_Rnd, db1Alpha, hFIG] = CBASS_U_ColoredNoiseMatchingRandSeries(db2LFP_Rnd, 1, inNOut, inNSampFit, true);
        CBASS_SaveFig(chOutPath, hFIG, {'ColoredNoiseFits'});
    end
    
    % Recomposes the signal and stores it in the output structure
    db2LFP_Rnd = (db2LFP_Rnd/db2Coeff)';
    
    % Check wether the size is good
    if inNOut ~= size(db2LFP_Rnd, 2), fprintf('inNOut = %d\t size(db2LPF_Rnd, 2) = %d\r', inNOut, size(db2LFP_Rnd, 2));
    else, fprintf('All good\r'); end
    
    % Plots the figure
    hFIG = figure('Position', [50 50 1250 600]);
    chFigName = ['Power_And_Example_', cMETHOD{iMth}];
    hPLT(1) = subplot(2, 2, 1);  CBASS_Plot_LFP_FourierPower(db2LFP_In, sREC.inSampleRate)
    title('LFP');
    hPLT(2) = subplot(2, 2, 2);  CBASS_Plot_LFP_FourierPower(db2LFP_Rnd, sREC.inSampleRate)
    title(sprintf('%s Surrogate LFP', cMETHOD{iMth}));
    linkaxes(hPLT, 'y');
    hPLT(1) = subplot(2, 2, 3); CBASS_Plot_LFP_EventExample(db2LFP_In, 1, [0 4], sREC.inSampleRate);
    hPLT(2) = subplot(2, 2, 4); CBASS_Plot_LFP_EventExample(db2LFP_Rnd, 5, [0 4], sREC.inSampleRate);
    linkaxes(hPLT, 'y'); pause(.1)
    CBASS_SaveFig(chOutPath, hFIG, {chFigName});
    clear hPLT
    
    % Checks the scaling
    db1Bin      = floor(min([db2LFP_In(:); db2LFP_Rnd(:)])*10)/10:0.1:ceil(max([db2LFP_In(:); db2LFP_Rnd(:)])*10)/10;
    db1H_In     = hist(db2LFP_In(:), db1Bin)./numel(db2LFP_In);
    db1H_Rnd    = hist(db2LFP_Rnd(:), db1Bin)./numel(db2LFP_Rnd);
    hFIG = figure; bar(db1Bin, [db1H_In' db1H_Rnd']);
    chFigName = ['Value_Histogram_' cMETHOD{iMth}];
    legend('Real', 'Surrogate'); title('LFP values');
    CBASS_SaveFig(chOutPath, hFIG, {chFigName});
end
end