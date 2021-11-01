function [db1Phase, bl1PulseCycle] = CBASS_U_PulsePhase(db2LFP, inSampleRate, inRefChan, db1Band, in1PulseIdx)
%CBASS Utility. Calculates a phase for a set of pulse. This is useful to estimate if points processes
%such as neuronal firing are phase locked to pulse using metrics like the PairWise Phase Concistency.
%Synopsis:
%	[db1Phase, bl1PulseCycle] = CBASS_U_PulsePhase(db2LFP, inSampleRate, inRefChan, db1Band, in1PulseIdx)
%
%Inputs:
%------------------------------------------------------------------------------------------------------------
%-db2LFP: 		A multidimensional (channel x time sample) array representing the signal of interest
%-inSampleRate: A positive number representing the sample rate of db2LFP.
%-inRefChan:	An integer indexing the channel of reference in db2LFP. This channel will be used to compute
% 				the phase.
%-db1Band: 		A two element row vector indicating the lower and upper bound of the band used to filter the 
% 				signal (example: [30 80]).
%-in1PulseIdx: 	An array of integer indexing the time samples at wich pulses of activity are observed. The
%				pulses are the output of CBASS.
%
%Outputs:
%------------------------------------------------------------------------------------------------------------
%-db1Phase: 		A row vector having as many elements as times samples in db2LFP. The vector gives the 
% 					phase of the signal in the reference channel db2LFP(inRefChan, :) when band pass filtered
%					using db1Band.
%-bl1PulseCycle: 	A row boolean vector having as many elements as times samples in db2LFP. The vector 
% 					indexes cycles around a pulse. This allows to select points in or outside a pulse 
% 					compare their phase.

%Check input:
if length(size(db2LFP)) > 2, error('db2LFP must be a matrix'); end
[inNChan, inNSmp] = size(db2LFP);
if any(size(inSampleRate) ~= 1)
	error('inSampleRate must be a single positive number') 
end
if inSampleRate < 0 
	error('inSampleRate must be a single positive number')
end
if any(size(inRefChan) ~= 1)
	error('inRefChan must index a single row in db2LFP') 
end
if mod(inRefChan, 1) ~= 0 || inRefChan < 1 || inRefChan > inNChan
	error('inRefChan must index a single row in db2LFP')
end
if length(db1Band) > 2 || any(diff(db1Band) <= 0) || any(db1Band < 0)
	error('db1Band must be a two element positive vector indexing a frequency band');
end
if any(mod(in1PulseIdx, 1) ~= 0) || any(in1PulseIdx < 0) || any(in1PulseIdx > inNSmp)
	error('in1PulseIdx must index columns in db2LFP');
end 

%Filters the LFP
[B, A] = butter(2, 2 * db1Band / inSampleRate);
db2_Filt_LFP = filtfilt(B, A, db2LFP(inRefChan, :)); %The LFP is transposed because filtfilt and hilbert opperates over rows

%Uses the Hilbert Transform
db2Hilbert 	= hilbert(db2_Filt_LFP);

%OUTPUT 1: Computes the phase
db1Phase   	= angle(db2Hilbert);

%Computes a vector indexing the cycles of each phase
db1Cycle 	= ceil(unwrap(db1Phase) ./  (2 * pi));

%Computes the cycles to which each pulse belongs
db1PulseCycle 	= db1Cycle(in1PulseIdx);

%OUTPUT 2: Computes a boolean indexing all cycles containing a pulse
bl1PulseCycle 	= ismember(db1Cycle, db1PulseCycle); 
