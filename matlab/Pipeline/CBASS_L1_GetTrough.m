function sTROUGH = CBASS_L1_GetTrough(db2LFP, inSampleRate, db1FilterBand, inRefChan, chLabel, chDataFormat)
% L1 of the bout pipeline: Identifies events in the activity band defined
% by db1FilterBand. Events correspond to the troughs of oscillatory
% activity at the band of interest in a reference channel inRefChan. They
% are represented as the Hilbert transform of the band filtered activity of
% each channel. The representation can either use complex or polar
% coordinates depending on the value of the optional variable chDataFormat.
%
% Input -------------------------------------------------------------------
%
% db2LFP:           a (channel x time sample) matrix containing the signal
%                   of interest - originaly meant to be a recording of the
%                   Local Field Potential(LFP) but it can be any
%                   multichannel time series.
% inSampleRate:     a positive number describing the sample rate of the
%                   signal
% db1FilterBand:    a (1 x 2) array describing the frequency band of
%                   interest i.e. [30 80] for 30 to 80 Hz.
% inRefChan:        a number specifying a reference channel. Events will be
%                   aligned to the trought of the band specific activity in
%                   that channel for trough identification. Default is the
%                   last channel of db2LFP.
% chLabel           (optional) a charactaer array describing the band of
%                   interest - (i.e. 'gamma' for [30 80Hz])
% chDataFormat:     (optional) a character array specifying the format of
%                   the hilbert transforms output. Can be 'complex' or
%                   'polar'. Default is 'complex'.
%
% Output ------------------------------------------------------------------
%
% sTROUGH:      a structure containing the following fields:
%               -.db1FilterBand an (1 x 2) array describing the frequency
%               band of interest i.e. [30 80] for 30 to 80 Hz.
%               -.db2Trough  a (2 * channel x trough) matrix containing the
%               hilbert transform of each channel of sREC.db2LFP filtered
%               in the band defined in sTROUGH.db1FilterBand at the trough
%               of the filtered signal in the reference channel inRefChan
%               -.in1Index the indices of the trough in sREC.db2LFP

% checks for the proper number of arguments
narginchk(4, 6);

% Creates a Label if not provided
if ~exist('chLabel', 'var'), chLabel = sprintf('%d-%dHz', db1FilterBand(1), db1FilterBand(2)); end
if ~exist('chDataFormat', 'var'), chDataFormat = 'complex'; end
if ~ismember(chDataFormat, {'polar', 'complex'})
    fprintf('%s is not a valid method set to umap', chDataFormat);
    chDataFormat = 'complex';
end

% Filters the LFP
[B, A] = butter(2, 2 * db1FilterBand / inSampleRate);
db2_Filt_LFP = filtfilt(B, A, db2LFP'); %The LFP is transposed because filtfilt and hilbert opperates over rows

% Uses the Hilbert Transform
db2_Hilbert = hilbert(db2_Filt_LFP);

% Computes the amplitude and the phase
db1_Amp     = abs(db2_Hilbert);
db1_Phase   = angle(db2_Hilbert);

% Finds the indices of troughs on the reference channels
bl1RefTrough = db1_Phase(1:end-1, inRefChan) > 0 & db1_Phase(2:end, inRefChan) < 0;

if strcmp(chDataFormat, 'complex')
    % Format the troughs so that each row is a motif and col 1:15 is the
    % real part and col 16:30 the imaginary part of the hilbert transform
    % of the filtered signal
    db2Trough       = [real(db2_Hilbert(bl1RefTrough, :)) imag(db2_Hilbert(bl1RefTrough, :))];
else
    % Formats the troughs so that that each row is a motif (the ref channel
    % is at the trough of the oscillation) and that col 1:15 correspond to
    % the amplitude of each channel and col 16:30 correspond to the phase.
    % Note these two groups should be normalized separately
    db2Trough       = [db1_Amp(bl1RefTrough, :) db1_Phase(bl1RefTrough, :)];
end

% Store the output
sTROUGH.chLabel         = chLabel;
sTROUGH.db1FilterBand   = db1FilterBand;
sTROUGH.db2Trough       = db2Trough;
sTROUGH.in1Index        = find(bl1RefTrough);
sTROUGH.chDataFormat    = chDataFormat;