function db2Filter = CBASS_U_MakeFilter(db2LFP, inSampleRate, in1EventIdx, db1Filter_WinSec)
% L2 of the bout pipeline: Identifies pulses of activity using the average
% of the LFP around a selctions of the troughs

% Gets the size of the LFP
[inNChan, inNSample] = size(db2LFP);

% Defines the window
in1Filt_RelIdx      = round(inSampleRate * db1Filter_WinSec(1)) : ...
    round(inSampleRate .* db1Filter_WinSec(2));

% Removes events if they are too close to the edge
bl1Rem              = in1EventIdx <= -in1Filt_RelIdx(1) | in1EventIdx >= inNSample - in1Filt_RelIdx(end);
in1EventIdx(bl1Rem) = [];
inNEvt              = length(in1EventIdx);

% Creates the ETA
db3LFP_ETA = nan(inNChan, length(in1Filt_RelIdx), inNEvt);
for iEvt = 1:inNEvt
    db3LFP_ETA(:, :, iEvt) = db2LFP(:, in1EventIdx(iEvt) + in1Filt_RelIdx);
end
db2Filter = nanmean(db3LFP_ETA, 3);