function db2LFP = PCP_U_NormalizeLFP(db2LFP, inSampleRate, varargin)
% Utility meant to normalize the LFP and make sure that it is independant
% of the inpendance of contacts

% Handles variable argument
if nargin > 2, in1ChanSel = varargin{1}; else, in1ChanSel = 1:size(db2LFP, 1); end

% Selects the channels
db2LFP = db2LFP(in1ChanSel, :);

% Defines a high pass filter
dbHighBound = 1;
[B, A] = butter(2, 2*dbHighBound/inSampleRate, 'high'); 

% Filters the LFP
db2LFP = filtfilt(B, A, db2LFP')';

% ZScores the LFP;
dbMu    = nanmean(db2LFP(:));
dbSig   = nanstd(db2LFP(:));
db2LFP = (db2LFP - dbMu)./dbSig;