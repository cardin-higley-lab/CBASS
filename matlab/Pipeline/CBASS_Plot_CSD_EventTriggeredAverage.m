function db2CSD_ETA_Inter = CBASS_Plot_CSD_EventTriggeredAverage(db2LFP, in1Event, db1WinSec, inSampleRate)
%Plotting utility of the bout pipeline to plot CSD event triggered averages

% Gets the number of channels and samples of the LFP
[inNChan, inNSample] = size(db2LFP);

% Makes sure that the event is in the index format
if islogical(in1Event), in1Event = find(in1Event); end
% Checks that the event vector has the right length
if any(mod(in1Event, 1) ~= 0) || any(in1Event) < 1 || any(in1Event) > inNSample
    error('in1Event must a vector of positive indices')
end

% Computes the relative indicies of the ETA
in1ETA_RelIdx   = round(inSampleRate * db1WinSec(1)): round(inSampleRate .* db1WinSec(2));

% Removes events if they are too close to the edge
bl1Rem          = in1Event <= -in1ETA_RelIdx(1) | in1Event >= inNSample - in1ETA_RelIdx(end);
in1Event(bl1Rem) = [];
inNEvt          = length(in1Event);

% Creates the ETA
db3LFP_ETA = nan(inNChan, length(in1ETA_RelIdx), inNEvt);
for iEvt = 1:inNEvt
    db3LFP_ETA(:, :, iEvt) = db2LFP(:, in1Event(iEvt) + in1ETA_RelIdx);
end
db2LFP_ETA = nanmean(db3LFP_ETA, 3);

% Sets the interpolation for the CSD
inInterp    = 20;
in1ItrIdx   = interp1(2:inNChan - 1, 2:inNChan - 1, 2:1/inInterp:inNChan - 1);

% Computes the CSD of the ETA
db2CSD_ETA = -diff(db2LFP_ETA, 2, 1);
db2CSD_ETA_Inter = interp1(2:inNChan - 1, db2CSD_ETA, in1ItrIdx);

% Plot the CSD
imagesc(in1ETA_RelIdx ./ inSampleRate, in1ItrIdx, db2CSD_ETA_Inter);
xlabel('Time (s)'); ylabel('Channel');