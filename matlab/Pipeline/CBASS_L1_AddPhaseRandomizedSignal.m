function sREC = CBASS_L1_AddPhaseRandomizedSignal(sREC, inSeed)
% Part of L1 of the bout pipeline. Creates surrogate LFP having the same
% covariance matrix and the same spectral density as sREC.db2LFP.
%
% Input -------------------------------------------------------------------
%
% sREC:         a structure requiring the following fields:
%               -.db2LFP a (channel x time sample) matrix containing the
%               signal (i.e. time series) of interest. 
%               -.inSampleRate a positive number representing the sample
%               rate of the time series.
% inSeed        (optional) a non negative integer used to seed the random
%               number generator.
% 
% Output ------------------------------------------------------------------
%
% sREC          the input structure with the additonal subfield:
%               -.db2LFP_Rnd a (channel x time sample) matrix of the same
%               size as sREC.db2LFP containing a phase randomized signal
%               having the same spectral density and the same covariance
%               matrix as sREC.db2LFP.

% Checks for the optional arguments
narginchk(1, 2);
if nargin < 2, inSeed = 1949; end

% Initialize the random number generator
rng(inSeed);

% Decomposes the LFP into orthogonal principal components
[db2Coeff, db2LFP_Rnd] = pca(sREC.db2LFP');

% Randomizes the phase of each principal component
db2LFP_Rnd = CBASS_U_PhaseRandomize1D(db2LFP_Rnd, 1);

% Recomposes the signal and stores it in the output structure
sREC.db2LFP_Rnd = (db2LFP_Rnd/db2Coeff)';