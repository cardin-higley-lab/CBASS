function chSigString = CBASS_Plot_Pulse_vs_State(bl1Pulse, bl1State, db1WinSec, inSampleRate, dbConvWinSec)
% Bout pipeline utility to estimate if the pulses in bl1Pulse are enriched
% during the state defined by bl1State. The function plots the average
% instantaneous frequency of pulse around the state ON transitions.

% Checks for optional argument
narginchk(4, 5);
if nargin < 5, dbConvWinSec = 0.5; end

% if bl1Pulse is not a logical vector assumes that it is a vector of
% indices. Converts it to a logical vector of the same size as bl1State
if ~islogical(bl1Pulse)
    in1Pulse = bl1Pulse;
    bl1Pulse = false(size(bl1State));
    bl1Pulse(in1Pulse) = true;
end

% Calculates the frequency of bouts and events during quiet and running 
dbFreq_StateON   = inSampleRate * sum(bl1Pulse & bl1State)./sum(bl1State);
dbFreq_StateOFF  = inSampleRate * sum(bl1Pulse & ~bl1State)./sum(~bl1State);

% Relative increase of frequency
dbRatio = dbFreq_StateON./dbFreq_StateOFF;
% Significance of the increase
dbPVal = 2 * (1 - normcdf(abs(dbFreq_StateON - dbFreq_StateOFF)./...
    sqrt((dbFreq_StateON ./ sum(bl1Pulse & bl1State)) +  (dbFreq_StateOFF ./ sum(bl1Pulse & ~bl1State)) )));

% Print the results on the screen
if dbPVal > 0.0001, chP = sprintf('%.4f', dbPVal); else, chP = '<0.0001'; end
chSigString = sprintf('OFF: %.2fHz\tON: %.2fHz\tRatio:%.2f\t(p = %s)\r', ....
    dbFreq_StateOFF, dbFreq_StateON, dbRatio, chP);

% Convolutes the trace by a square window to compute an instantaneous frequency
db1BoutTrace    = conv(bl1Pulse, rectwin(dbConvWinSec * inSampleRate), 'same')./dbConvWinSec;

% Defines the relative indices of the ETA 
in1ETA_RelIdx   = round(inSampleRate * db1WinSec(1)): round(inSampleRate .* db1WinSec(2));

% Sets the event set
in1RunON            = 1 + find(~bl1State(1:end-1) & bl1State(2:end));
bl1Rem              = in1RunON <= -in1ETA_RelIdx(1) | in1RunON >= length(bl1State) - in1ETA_RelIdx(end);
in1RunON(bl1Rem)    = [];

%  Calculates the event triggered avergae
db1ETA_Bout     = mean(db1BoutTrace(in1ETA_RelIdx + in1RunON'));
db1ETSEM_Bout   = std(db1BoutTrace(in1ETA_RelIdx + in1RunON'))./sqrt(length(in1RunON));

% Defines the time
db1Time     = in1ETA_RelIdx./inSampleRate;

% Plot the results
fill([db1Time db1Time(end:-1:1)], ...
    [db1ETA_Bout + db1ETSEM_Bout  db1ETA_Bout(end:-1:1) - db1ETSEM_Bout(end:-1:1)], ...
    [.5 .5 .5], 'LineStyle', 'none', 'FaceAlpha', .3); hold on
plot(db1Time, db1ETA_Bout, 'k'); hold on
% plot(db1Time, (db1ETA_Bout + [db1ETSEM_Bout; -db1ETSEM_Bout])', 'k--');
db1YL = ylim; plot([0 0], db1YL, 'r--');
ylabel('Instantaneous Frequency (Hz)'); xlabel('Time (s)');