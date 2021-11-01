function [sFREQ_BAND, cTROUGH, cTRGH_RND, cFILTER, cPULSE] = ...
    CBASS_Main_DetectEnrichedMotif(db2LFP, inSampleRate, cBAND, cSTATE, sOPTION)
% Deprecated main function for the detection of enriched band specific
% activity motif. Chance levels of detection and significance levels are
% estimated by repeating procedures on surrogate data having the same
% spectral density the same covariance matrix between channels. The
% surrogate data are generated with the function:
% CBASS_L1_AddPhaseRandomizedSignal. The function makes use of 3 additional
% subroutines:
%
% CBASS_L1_GetTrough:    performs trough identification of real and
%                       surrogate signals using the Hilbert transform.
% CBASS_L2_GetFilters:   dentifies groups of troughs enriched for the
%                       state of interst. Uses the average LFP activity 
%                       around theme to generate template motifs 
%                       (i.e. filters). And clusters filters based on their
%                       similarity (i.e. correlation) using spectral
%                       clustering.
% CBASS_L3_GetPulse_2:   find events in the signal (i.e. pulses) matching
%                       the template motifs
% 
%
% Input -------------------------------------------------------------------
% db2LFP:       a (channel x time sample) matrix containing the signal of
%               interest - originaly meant to be a recording of the Local 
%               Field Potential(LFP) but it can be any multichannel time 
%               series.
% inSampleRate: a positive number describing the sampling rate of the
%               signal
% cBAND:        a (1 x 2) array describing the frequency band of interest
%               i.e. [30 80] for 30 to 80 Hz OR several such arrays
%               grouped in a cell array. The analysis will be performed
%               for each element.
% cSTATE:       a logical vector, containing as many elements as time
%               samples in db2LFP, indexing the state in which enriched
%               band specific activity is observed OR a cell array of such
%               vector. If so the number of element of the cell array must
%               mach the number of element in cBAND.
%
% Option structure --------------------------------------------------------
% The variable sOPTION can be used to pass optional arguments. It may
% contain the following fields:
%
% .cBAND_LABEL:     a cell array of size matching cBAND containings labels
%                   for the band of interest - (i.e. 'gamma' for [30 80Hz])
% .chDataFormat:    a character array specifying the format of the hilbert
%                   transforms output. Can be 'complex' or 'polar'. Default
%                   is 'complex'.
% .inRefChan:       a number specifying a reference channel. Events will be
%                   aligned to the trought of the band specific activity in
%                   that channel for trough identification. Default is the
%                   last channel of db2LFP.
% .cBASELINE        a logical vector, containing as many elements as time
%                   samples in db2LFP, indexing the state in which enriched
%                   band specific activity is NOT observed OR a cell array
%                   of such vector. If so the number of element of the cell
%                   array must mach the number of element in cBAND.
% .blZScore:        logical specifying if trought data is to be zscored for 
%                   k-means partitioning. Default is true.
% .inMethod:        single integer with value 1 or 2. Determines the method
%                   used to set the threshold of similarity used to retain
%                   an edge in spectral clustering during the grouping
%                   of templates. Method 1 sets a general threshold as
%                   the average similarity of templates identified on
%                   surrogate data. Method 2 sets a threshold for each edge
%                   as the maximum simarity between its nodes and a chance
%                   template computed as the average activity around all
%                   troughs in the surrogate data. Method 2 identifies
%                   more templates, but these templates will generate
%                   overlapping sets of pulses. Default is 1.
% .inNClu:          number of cluster used for k-means partitioning.
%                   Default is 20
% .dbSigThrs:       threshold for the significance of the enrichment in 
%                   trough partition.  P-Values are computed with a 
%                   binomial test. Default is 10.^-4.
% .inNMaxIter:      maximum iteration used by the k-means algorithm. 
%                   Default is 10000.
% .blVerbose:       logical setting whether processing updates should be
%                   displayed in the command window. Default is true.
% .blCntr:          logical setting whether data and template are centered 
%                   about their mean in template matching. Default is true.
% .blNorm:          logical setting whether data should be normalized by 
%                   its sliding S.D. during template matching. Default is 
%                   true.
%
% Output ------------------------------------------------------------------
% sFREQ_BAND:       a structure array of the same size as cBAND containing
%                   the following fields:
%   .db1Band        The value of cBAND for that instance of the array
%   .chBandLabel    The label of the band for that instance of the array 
%                   (may be specified in sOPTION.cBAND_LABEL).
%   .in1TroughIdx   the indices of the trough of oscillatory activity in
%                   sREC.db2LFP used to define activity templates (Troughs
%                   are generated by CBASS_L1_GetTrough).                  
%   .sMOTIF         a structure array whose length equals the number of
%                   template identified and containing the following
%                   fields:
%                   -.bl1TroughSel a logical vector indexing the trough
%                   used to compute the template motif
%                   -.db2Filter a (channel x time sample) matrix describing
%                   the template motif. The number of time samples is one
%                   period of the median frequency of the band of interest.
%                   -.db1TM_Score the score of the template matching. A
%                   vector having as many elements as time samples in
%                   db2LFP.
%                   -.dbSigTheshold the theshold for significance for
%                   peaks in .db1TM_Score (see CBASS_L3_GetPulse_2 for
%                   details)
%                   -.bl1Pulse a boolean, having as many elements as time
%                   samples in db2LFP, indexing significant pulses of band
%                   specific activity. <---  FINAL OUTPUT OF THE PROCEDURE.
%
% Optional Output ---------------------------------------------------------
% cTROUGH           a cell array containing the output of the subroutine
%                   CBASS_L1_GetTrough for each frequency band in cBAND.
%                   Each instance of cTROUGH is a structure sTROUGH
%                   containing the following fields:
%   .db1FilterBand  an (1 x 2) array describing the frequency
%                   band of interest i.e. [30 80] for 30 to 80 Hz.
%   .db2Trough      a (2 * channel x trough) matrix containing the
%                   hilbert transform of each channel of sREC.db2LFP filtered
%                   in the band defined in sTROUGH.db1FilterBand at the trough
%                   of the filtered signal in the reference channel inRefChan
%   .in1Index       the indices of the trough in sREC.db2LFP
%
% cTRGH_RND         a cell array containing the output of the subroutine
%                   CBASS_L1_GetTrough for each frequency band in cBAND
%                   applied to the surrogate signal generated by
%                   CBASS_L1_AddPhaseRandomizedSignal. Each instance of
%                   cTRGH_RND is a structure sTRGH_RND containing the same
%                   fields as the elements of cTROUGH (see above)
%
% cFILTER           a cell array containing the output of the subroutine
%                   CBASS_L2_GetFilters for each frequency band in cBAND.
%                   Each instance of cFILTER is a structure array sFILTER
%                   where each element is a distinct filter (i.e. template)
%                   containing the following fields:
%   .in1CluSel      a vector containing the indices of the k-means
%                   clusters used to build the filter
%   .bl1Member      a boolean indexing troughs used to build the
%                   template motif
%   .inNObs         the number of troughs in sFILTER.bl1Member.
%   .dbRate         the fraction of the troughs occuring during the
%                   state of interest.
%   .dbRate_Dev     the enrichment for the state of interest i.e.
%                   sFILTER.dbRate minus the overall ratio of troughs
%                   occuring during the state of interest
%   .dbPVal         the p-value of a binomial test of the of dbRate
%                   compared to the overall rate of occurence of the state
%                   of interest.
%   .db2Filter      a (channel x time sample) matrix describing the
%                   template motif. The number of time samples is one
%                   period of the median frequency of the band of interest.
%
% cPULSE            a cell array containing the output of the subroutine
%                   CBASS_L3_GetPulse_2 for each frequency band in cBAND.
%                   Each instance of cPULSE is a structure sPULSE where
%                   each element is generated with a distinct filter (i.e.
%                   template) and containing the following fields:
%   .db1Score       the score of the template matching. A vector having
%                   as many elements as time samples in sREC.db2LFP.
%   .bl1Peak        a boolean indexing the local maxima in db1Score
%   .db1Score_H0    the score obtained on surrogate data
%   .bl1Peak_H0     a boolean indexing local maxima on db1Score_H0
%   .dbThrPeak      the theshold for significance for peaks.
%                   Significant peaks are considered pulses of band
%                   specific activity. Set as the 95% percentile of
%                   sPULSE.db1Score_H0(sPULSE.bl1Peak_H0)
%   .bl1Pulse       a boolean indexing significant pulses of band
%                   specific activity. It has as many elements as time
%                   samples in sREC.db2LFP.     <----- FINAL OUTPUT OF THE PROCEDURE.
%   .dbP_KS         the p-value of a Kolmogorov-Smirnov test of the
%                   difference between the distrubutions of score values at
%                   peak in the LFP and in the surrogate data
%   .dbP_Rate       the p-value of a binomial test of the occurence of
%                   significant pulses in the real data compared to the
%                   surrogate data (where it is set as 5% by definition -
%                   see field .db1ThrPeak).

%Checks the number of arguments -------------------------------------------
narginchk(4, 5);
if nargin < 5, sOPTION = struct; end

%Checks for non optional arguments ----------------------------------------
[cBAND, cSTATE, inNChan, bl1Remove] = CheckArg(db2LFP, inSampleRate, cBAND, cSTATE); % Function at the end of the script

% Deals with the option structure -----------------------------------------
sOPTION = CheckOption(sOPTION, cBAND, cSTATE, inNChan, bl1Remove); % Function at the end for the script

% Computes ----------------------------------------------------------------
    %Formats the LFP and computes the phase randomized signal (this will be
    %used for chance level estimation and statistical testing)
if sOPTION.blVerbose, fprintf('---->> Compute phase randomized signal ... '); end
sREC.db2LFP         = db2LFP;
sREC.inSampleRate   = inSampleRate;
sREC                = CBASS_L1_AddPhaseRandomizedSignal(sREC);
if sOPTION.blVerbose, fprintf('Done\r'); end

    %Initializes the output structure
inNBnd = length(cBAND); 
sFREQ_BAND = struct('db1Band', cell(1, inNBnd), 'chBandLabel', cell(1, inNBnd), ...
    'in1TroughIdx', cell(1, inNBnd), 'sMOTIF', cell(1, inNBnd));

    %Initializes the optional output cell arrays
[cTROUGH, cTRGH_RND, cFILTER, cPULSE] = deal(cell(size(cBAND))); 

    %Loops through bands of interest 
for iBnd = 1:inNBnd
    % Choose the state of interest
    if sOPTION.blVerbose, fprintf('\r------ %s ---------------------------\r', sOPTION.cBAND_LABEL{iBnd}); end
    bl1State    = cSTATE{iBnd};

    % Extracts trougths for the real and surrogates signals
    if sOPTION.blVerbose, fprintf('---->> Extract hilbert troughs ... '); end
    sTROUGH     = CBASS_L1_GetTrough(sREC.db2LFP, inSampleRate, cBAND{iBnd}, ...
        sOPTION.inRefChan, sOPTION.cBAND_LABEL{iBnd}, sOPTION.chDataFormat);
    sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, inSampleRate, cBAND{iBnd}, ... 
        sOPTION.inRefChan, sOPTION.cBAND_LABEL{iBnd}, sOPTION.chDataFormat);
    if sOPTION.blVerbose, fprintf('Done\r\r'); end
    
    % Gets the filters
    if sOPTION.blVerbose, fprintf('---->> Compute filters ... \r'); end
    sFILTER =   CBASS_L2_GetFitlters(sREC, sTROUGH, sTRGH_RND, bl1State, ...
        sOPTION.cBASELINE{iBnd}, sOPTION.blZScore, sOPTION.inMethod, sOPTION.inNClu, sOPTION.dbSigThrs, ...
        sOPTION.inNMaxIter, sOPTION.blVerbose);
    
    % Get pulses for each filter
    if sOPTION.blVerbose, fprintf('\r---->> Detect pulses ... '); end
    inNFlt  = length(sFILTER); 
    sMOTIF = struct('bl1TroughSel', cell(1, inNFlt), 'db2Filter', cell(1, inNFlt), ...
        'db1TM_Score', cell(1, inNFlt), 'bl1Pulse', cell(1, inNFlt)); % Intialize the final motif structure
    sPULSE = struct('db1Score', cell(1, inNFlt), 'bl1Peak', cell(1, inNFlt), ...
        'db1Score_H0', cell(1, inNFlt), 'bl1Peak_H0', cell(1, inNFlt), 'dbThrPeak', cell(1, inNFlt), ...
        'bl1Pulse', cell(1, inNFlt), 'dbP_KS', cell(1, inNFlt), ...
        'dbP_Rate', cell(1, inNFlt)); % Initializes the L3 pulse structure
    
    % Loops through filters
    for iFlt = 1:inNFlt
        sPULSE(iFlt) = CBASS_L3_GetPulse_2(sREC, sFILTER(iFlt).db2Filter, ...
            sOPTION.blCntr, sOPTION.blNorm);
        % Fills the output structure
        sMOTIF(iFlt).bl1TroughSel   = sFILTER(iFlt).bl1Member; % Subset of the troughs used to compute the template
        sMOTIF(iFlt).db2Filter      = sFILTER(iFlt).db2Filter; % Template used of template matching (i.e. filter)
        sMOTIF(iFlt).db1TM_Score    = sPULSE(iFlt).db1Score; % Score for template matching
        sMOTIF(iFlt).bl1Pulse       = sPULSE(iFlt).bl1Pulse; % logical indexing pulses in the input signal
    end
    if sOPTION.blVerbose, fprintf('Done\r\r'); end
    
    % Aggregates pulse data
    sFREQ_BAND(iBnd).db1Band        = cBAND{iBnd};
    sFREQ_BAND(iBnd).chBandLabel    = sOPTION.cBAND_LABEL{iBnd};
    sFREQ_BAND(iBnd).in1TroughIdx   = sTROUGH.in1Index;
    sFREQ_BAND(iBnd).sMOTIF         = sMOTIF;
    
    % Aggregate intermediary steps
    cTROUGH{iBnd}   = sTROUGH;
    cTRGH_RND{iBnd} = sTRGH_RND;
    cFILTER{iBnd}   = sFILTER; 
    cPULSE{iBnd}    = sPULSE;
end

%------------------------------------------------------------------------------------------------  
function [cBAND, cSTATE, inNChan, bl1Remove] = CheckArg(db2LFP, inSampleRate, cBAND, cSTATE)
%Utility to check non optional arguments

% Checks that the signal is a 2D matrix
if ~ismatrix(db2LFP); error('The signal must be a (channel x time sample) matrix'); end
[inNChan, inNSamp] = size(db2LFP);

% Check sample rate
if ~(numel(inSampleRate) > 1 || any(mod(inSampleRate(:), 0) ~= 0))
    error('The sample rate must be a single positive number')
end

% Checks if states and band are not cell array
if ~iscell(cBAND); cBAND = {cBAND}; end
if ~iscell(cSTATE); cSTATE = {cSTATE}; end

% Checks that cBAND and cSTATE are properly formatted and have the same size
if ~isvector(cBAND); error('cBAND must be a vector cell array'); end
if ~isvector(cSTATE); error('cSTATE must be a vector cell array'); end
if iscolumn(cBAND); cBAND = cBAND'; end
if iscolumn(cSTATE); cSTATE = cSTATE'; end
if any(size(cBAND) ~= size(cSTATE)); error('cSTATE and cBAND must have the same size'); end

% Checks if on the entry of the state and bands is properly formated
bl1BandCorrect = cellfun(@(x) all(size(x) == [1 2]) & all(x > 0), cBAND);
if any(~bl1BandCorrect)
    fprintf('Some instances of cBAND are improperly formatted and will be skipped\r');
end
bl1StateCorrect = cellfun(@(x) isvector(x) & islogical(x) & numel(x) == inNSamp, cSTATE);
if any(~bl1StateCorrect)
    fprintf('Some instances of cSTATE are improperly formatted and will be skipped\r');
end
    
% Skips improperly formatted input
bl1Remove   = ~bl1BandCorrect | ~bl1StateCorrect;
cBAND(bl1Remove)    = [];
cSTATE(bl1Remove)   = [];
if isempty(cSTATE); error('No valid entry for state and/or frequency bands'); end

%------------------------------------------------------------------------------------------------
function sOPTION = CheckOption(sOPTION, cBAND, cSTATE, inNChan, bl1Remove)
% Utility to check the options structure. Makes use of the utility
% CheckField present at the end of the script

% Checks the band label
blLblErr = false;
if ~isfield(sOPTION, 'cBAND_LABEL')
    sOPTION.cBAND_LABEL = cellfun(@(x) sprintf('%d-%dHz', x(1), x(2)), cBAND, 'UniformOutput' , false); 
elseif ~iscell(sOPTION.cBAND_LABEL)
    if ischar(sOPTION.cBAND_LABEL), sOPTION.cBAND_LABEL = {sOPTION.cBAND_LABEL};
    else, sOPTION.cBAND_LABEL = cellfun(@(x) sprintf('%d-%dHz', x(1), x(2)), cBAND, 'UniformOutput' , false); 
        blLblErr = true; 
    end
elseif ~all(size(sOPTION.cBAND_LABEL) == size(bl1Remove))   
    sOPTION.cBAND_LABEL = cellfun(@(x) sprintf('%d-%dHz', x(1), x(2)), cBAND, 'UniformOutput' , false); blLblErr = true;
else
    sOPTION.cBAND_LABEL(bl1Remove) = [];
    in1Wrong = find(~cellfun(@ischar, sOPTION.cBAND_LABEL));
    for iWrg = in1Wrong
        sOPTION.cBAND_LABEL{iWrg} = sprintf('%d-%dHz', cBAND{iWrg}(1), cBAND{iWrg}(2));
        fprintf('Instance %d of sOPTION.cBAND_LABEL was set to default\r', iWrg);
    end
end
if blLblErr, fprintf('sOPTION.cBAND_LABEL is not valid. Set to default\r'); end

% Checks L1 options for formatting hilbert troughs
sOPTION = CheckField(sOPTION, 'chDataFormat', @(x) ~ismember(x, {'polar', 'complex'}), 'complex');
sOPTION = CheckField(sOPTION, 'inRefChan', @(x) numel(x) > 1 || any(mod(x(:), 1) ~= 0) || any(x(:) > inNChan), inNChan);

% Checks L2 options for the computation of filters
    % Checks baseline option argument
blBasErr = false;
if ~isfield(sOPTION, 'cBASELINE')
    sOPTION.cBASELINE = cellfun(@(x) ~x, cSTATE, 'UniformOutput' , false);
elseif ~iscell(sOPTION.cBASELINE)
    if islogical(sOPTION.cBASELINE) & all(size(sOPTION.cBASELINE) == cSTATE{1}), sOPTION.cBASELINE = {sOPTION.cBASELINE};
    else, sOPTION.cBASELINE = cellfun(@(x) ~x, cSTATE, 'UniformOutput' , false); blBasErr = true; end
elseif ~all(size(sOPTION.cBASELINE) == size(bl1Remove))
    sOPTION.cBASELINE = cellfun(@(x) ~x, cSTATE, 'UniformOutput' , false); blBasErr = true;
else
    sOPTION.cBASELINE(bl1Remove) = [];
    for iTst = 1:2
        if iTst == 1, in1Wrong = find(~cellfun(@(x, y) islogical(x) & all(size(x) == size(y)), sOPTION.cBASELINE, cSTATE));
        else, in1Wrong = find(cellfun(@(x, y) any(x & y), sOPTION.cBASELINE, cSTATE)); end
        for iWrg = in1Wrong
            sOPTION.cBASELINE{iWrg} = ~cSTATE{iWrg};
            fprintf('Instance %d of sOPTION.cBASELINE was set to default\r', iWrg);
        end
    end
end
if blBasErr, fprintf('sOPTION.cBASELINE does not match the format of cSTATE. Set to default\r'); end
    % Checks other optional arguments
sOPTION = CheckField(sOPTION, 'blZScore', @(x) numel(x) > 1 || any(~islogical(x(:))), true);
sOPTION = CheckField(sOPTION, 'inMethod', @(x) numel(x) > 1 || any(~ismember(x(:), [1 2])), 1);
sOPTION = CheckField(sOPTION, 'inNClu', @(x) numel(x) > 1 || any(mod(x(:), 1) ~= 0) || any(x(:) <= 0), 20);
sOPTION = CheckField(sOPTION, 'dbSigThrs', @(x) numel(x) > 1 || any(x(:) <= 0) || any(x(:) > 1), 10.^-4);
sOPTION = CheckField(sOPTION, 'inNMaxIter', @(x) numel(x) > 1 || any(mod(x(:), 1) ~= 0) || any(x(:) <= 0), 10000);
sOPTION = CheckField(sOPTION, 'blVerbose', @(x) numel(x) > 1 || any(~islogical(x(:))), true);

% Checks L3 options for template matching
sOPTION = CheckField(sOPTION, 'blCntr', @(x) numel(x) > 1 || any(~islogical(x(:))), true);
sOPTION = CheckField(sOPTION, 'blNorm', @(x) numel(x) > 1 || any(~islogical(x(:))), true);


function sOPTION = CheckField(sOPTION, chField, fCOND, DefaultVal)
% Utility to check the field of an option
if ~isfield(sOPTION, chField),  sOPTION.(chField) = DefaultVal; end
if fCOND(sOPTION.(chField))
    fprintf('sOPTION.%s is not valid. Set to default\r', chField);
    sOPTION.(chField) = DefaultVal;
end