function sREC = CBASS_L1_AddPhaseRandomizedSignal_2(sREC, bl1RefSel, blFit, inSeed)
% Alternate part of L1 of the bout pipeline. Creates surrogate LFP having
% the same covariance matrix and a spectral density matched to the LFP
% contained in sREC. The surrogate signal can be matched to a subselection
% of the reference LFP signal with the optional boolean bl1RefSel.
% bl1RefSel must be the same lenght as sREC.db2LFP. If the optional
% argument blFit is set to true a 1/f^alpha model is fitted to the
% reference data to generate the surrogate signal.
%
% Input -------------------------------------------------------------------
%
% sREC:         a structure requiring the following fields:
%               -.db2LFP a (channel x time sample) matrix containing the
%               signal (i.e. time series) of interest. 
%               -.inSampleRate a positive number representing the sample
%               rate of the time series.
% bl1RefSel     (optional) a logical vector having the same number of
%               elements as time samples in sREC.db2LFP indexing time
%               sample which will be used to derive the spectral density
%               that will be reproduced in the output surrogate signal.
%               Default: the full signal is used.
% blFit         (optional) if true the spectral density of the output
%               will follow a 1/f^alpha model fitted to the spectral
%               density of the reference signal. Default is false.
% inSeed        (optional) a non negative integer used to seed the random
%               number generator.
% 
% Output ------------------------------------------------------------------
%
% sREC          the input structure with the additonal subfield:
%               -.db2LFP_Rnd a (channel x time sample) matrix of the same
%               size as sREC.db2LFP containing a phase randomized signal
%               having a matched spectral density and the same covariance
%               matrix as sREC.db2LFP.

% Gets the number of samples
inNSmp = size(sREC.db2LFP, 2);

% Checks for the optional arguments
narginchk(1, 4);
if nargin < 2, bl1RefSel = true(1, inNSmp); end
if nargin < 3, blFit = false; end
if nargin < 4, inSeed = 1949; end

% Initialize the random number generator
rng(inSeed);

% Selects the reference LFP
db2LFP_Ref  = sREC.db2LFP(:, bl1RefSel);

% Decomposes the LFP into orthogonal principal components
[db2Coeff, db2LFP_Rnd] = pca(db2LFP_Ref');

% Randomizes the phase of each principal component
if blFit
    dbFreqMax   = 120;
    inNSmpFit   = round(dbFreqMax .* size(db2LFP_Ref, 2) ./ sREC.inSampleRate);
    db2LFP_Rnd  = CBASS_U_ColoredNoiseMatchingRandSeries(db2LFP_Rnd, 1, inNSmp, inNSmpFit);
else
    db2LFP_Rnd  = CBASS_U_PowerMatchingRandSeries(db2LFP_Rnd, 1, inNSmp);
end

% Recomposes the signal and stores it in the output structure
sREC.db2LFP_Rnd = (db2LFP_Rnd/db2Coeff)';