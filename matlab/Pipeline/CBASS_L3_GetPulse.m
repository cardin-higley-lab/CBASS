function sPULSE = CBASS_L3_GetPulse(sREC, sTROUGH, bl1Member, inFilt_NCycle, blCntr, blNorm)
% L3 of the bout pipeline: Identifies pulses of activity using a sliding
% cosine similarity between the signal contained in sREC.db2LFP and
% spatio-temporal motifs (i.e. filter). Filters are computed as the average
% signal activity around the subset of troughs of input structure sTROUGH
% indexed by the boolean bl1Member. Filters span a number of period of the
% average band used to define sTROUGH defined by the variable
% inFilt_NCycle. The best result will be obtained when the signal and
% template have been centered to their mean and normalized by their S.D..
% In the case the cosine similarity is equivalent to the correlation. The
% function computes the filteres using the average. 
%
% NOTE: This early version of the L3 computes template (i.e. filters) out
% of a subset of the events of the input structure sTROUGH indexed by the
% logical array bl1Member. It was initially meant to be used with the
% output of the deprecated function CBASS_L2_GetEnrichedRegion.m. It can
% also be used with the structure sCLU optionally returned by
% CBASS_L2_GetFilters. Finally, it can be used to test different lengths for
% the templates (expressed in number of cycles of the median frequency of
% the band of interst -- CBASS_L2_GetFilters returns template having one
% cycle in length -- results tend to degrade when increasing the number of
% cycles).
%
% Input -------------------------------------------------------------------
%
% sREC:         a structure requiring the following fields:
%               -.db2LFP a (channel x time sample) matrix containing the
%               signal (i.e. time series) of interest. 
%               -.db2LFP_Rnd a (channel x time sample) matrix containing
%               the surrogate signal (generated by
%               CBASS_L1_AddPhaseRandomizedSignal) --- Will be computed if
%               omitted.
%               -.inSampleRate a positive number representing the sample
%               rate of the time series.
% sTROUGH:      a structure requiring the following fields:
%               -.in1Index the indices of the trough in sREC.db2LFP
%               -.db1FilterBand the frequency band used to define the
%               troughs. Formatted as a 2 element row vector describing the
%               lower and upper bound of the band (i.e. [30 80] for 30 to
%               80Hz).
% inFilt_NCycle (optional) a positive number defining the number of periods
%               of the median frequency of the sTROUGH.db1FilterBand that
%               the filter will contain. Default is 1.
% blCntr        (optional) a boolean setting whether signal
%               and template should be centered about their mean. Default
%               is true.
% blNorm        (optional) a boolean setting whether signal and template
%               should be normalized by their S.D.. Default is true.
% 
% Output ------------------------------------------------------------------
%
% sPULSE:   a structure array containing the following fields:
%           -.db1Score the score of the template matching. A vector having
%           as many elements as time samples in sREC.db2LFP.
%           -.bl1Peak a boolean indexing the local maxima in db1Score
%           -.db1Score_H0 the score obtained on surrogate data
%           -.bl1Peak_H0 a boolean indexing local maxima on db1Score_H0
%           -.dbThrPeak the theshold for significance for peaks.
%           Significant peaks are considered pulses of band specific
%           activity. Set as the 95% percentile of
%           sPULSE.db1Score_H0(sPULSE.bl1Peak_H0)
%           -.bl1Pulse a boolean indexing significant pulses of band
%           specific activity. It has as many elements as time samples in
%           sREC.db2LFP. <----- FINAL OUTPUT OF THE PROCEDURE.
%           -dbP_KS the p-value of a Kolmogorov-Smirnov test of the
%           difference between the distrubutions of score values at peak in
%           the LFP and in the surrogate data 
%           -dbP_Rate the p-value of a binomial test of the occurence of
%           significant pulses in the real data compared to the surrogate
%           data (where it is set as 5% by definition - see field
%           .db1ThrPeak).

% Sets optional parameters if not provided
narginchk(3, 6);
if ~exist('inFilt_NCycle', 'var'), inFilt_NCycle = 1; elseif isempty(inFilt_NCycle), inFilt_NCycle = 1; end
if ~exist('blCntr', 'var'), blCntr = true; elseif isempty(blCntr), blCntr = true; end
if ~exist('blNorm', 'var'), blNorm = true; elseif isempty(blNorm), blNorm = true; end

% Gets the size of the LFP
[inNChan, inNSample] = size(sREC.db2LFP);

% Gets the indices of the epoches on the initial LFP recording
in1EventIdx     = sTROUGH.in1Index(bl1Member); % Clarifies event idx

% Defines the window
dbCycLen            = 1./mean(sTROUGH.db1FilterBand);
db1Filter_WinSec    = [-dbCycLen dbCycLen] * inFilt_NCycle / 2;
in1Filt_RelIdx      = round(sREC.inSampleRate * db1Filter_WinSec(1)) : ...
    round(sREC.inSampleRate .* db1Filter_WinSec(2));

% Removes events if they are too close to the edge
bl1Rem              = in1EventIdx <= -in1Filt_RelIdx(1) | in1EventIdx >= inNSample - in1Filt_RelIdx(end);
in1EventIdx(bl1Rem) = [];
inNEvt              = length(in1EventIdx);

% Creates the ETA
db3LFP_ETA = nan(inNChan, length(in1Filt_RelIdx), inNEvt);
for iEvt = 1:inNEvt
    db3LFP_ETA(:, :, iEvt) = sREC.db2LFP(:, in1EventIdx(iEvt) + in1Filt_RelIdx);
end
db2Filter = nanmean(db3LFP_ETA, 3);

% Runs template matching, detects the peaks and computes the enveloppe
db1Score    = CBASS_U_MultiChannelTemplateMatching(sREC.db2LFP, db2Filter, blCntr, blNorm);
bl1Peak     = [false, diff(db1Score(1:end - 1)) > 0 & diff(db1Score(2:end)) < 0, false];
% db1Env      = abs(hilbert(db1Score));

% Create a randomized LFP if not already present in the sREC structure
if ~isfield(sREC, 'db2LFP_Rnd')
    sREC = CBASS_L1_AddPhaseRandomizedSignal(sREC);
end

% Runs the template matching on the phase randomized LFP, detects the peaks
% and computes the enveloppe
db1Score_H0     = CBASS_U_MultiChannelTemplateMatching(sREC.db2LFP_Rnd, db2Filter, blCntr, blNorm);
bl1Peak_H0      = [false, diff(db1Score_H0(1:end - 1)) > 0 & diff(db1Score_H0(2:end)) < 0, false];
% db1Env_H0      = abs(hilbert(db1Score_Rnd));

% Detects significant peak: i.e. pulses
dbThrPeak   = quantile(db1Score_H0(bl1Peak_H0), .95);
bl1Pulse    = bl1Peak & db1Score > dbThrPeak;

% % Detects significant epochs i.e. bouts
% dbThrEnv    = quantile(db1Env_Rnd, .95);
% bl1Bout     = db1Env > dbThrEnv;

% Tests that the distribution is different and gives a warning if not
[blSig, dbP_KS] = kstest2(db1Score(bl1Peak), db1Score_H0(bl1Peak_H0));
if ~blSig, warning('Peak distribution is not different from noise'); end

% Test that the proportion of significant events is above chance
dbP_Rate = CBASS_U_TestRate(mean(db1Score(bl1Peak) > dbThrPeak), sum(bl1Peak), .05);
if dbP_Rate > .05; warning('Peak Rate not different from noise'); end

% Stores the result in the output structure
sPULSE.db1Score     = db1Score;
sPULSE.bl1Peak      = bl1Peak;
sPULSE.db1Score_H0  = db1Score_H0;
sPULSE.bl1Peak_H0   = bl1Peak_H0;
sPULSE.dbThrPeak    = dbThrPeak;
sPULSE.bl1Pulse     = bl1Pulse;
sPULSE.dbP_KS       = dbP_KS;
sPULSE.dbP_Rate     = dbP_Rate;