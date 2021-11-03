function [db2Power, in1CntrIdx] = CBASS_U_FourierPower(db1Trace, inSampleRate, varargin)
%[DB2POWER, IN1CNTRIDX] = CBASS_U_FourierPower(DB1TRACE, INSAMPLERATE,
%[DBWINLENSEC, INNSTEP, BLDOPLOT]) computes the spectral power of DB1TRACE
%over time. DB1TRACE is divided into segments whose length is specified by
%DBWINLENSEC and INSAMPLERATE. INNSTEP controls the overlap between
%segments (e.g. if INNSTEP is 4, then segments are spaced with 1/4 of the
%length of the window). Defautl for DBWINLENSEC is 0.5 sec. Default for
%INNSTEP is 4. The power spectrum is returned in DB2POWER. In addition, the
%indices of the center of each segment on DB1TRACE are returned in
%IN1CNTRIDX.
%
%The code is adapted from the help page 'power spectral density estimate
%using fft' and is the same as DoFourierPowerSpectrum
%
%2016-11-15 QP: Transcripted from DoFourierPowerSpectrum


switch length(varargin)
    case 0
        dbWinLenSec = 0.5;
        inNStep = 5;
        blDoPlot = false;
    case 1
        dbWinLenSec = varargin{1};
        inNStep = 5;
        blDoPlot = false;
    case 2
        dbWinLenSec = varargin{1};
        inNStep = varargin{2};
        blDoPlot = false;
    otherwise
        dbWinLenSec = varargin{1};
        inNStep = varargin{2};
        blDoPlot = varargin{3};
end

%Checks that step divides the traces without remainder
inWinLen    = round(dbWinLenSec * inSampleRate);
inStepLen   = inWinLen / inNStep;
if mod(inStepLen, 1) ~= 0
    error('The step taken does not divide the length into segments of equal length');
end

%Checks for the maximal number of chunks that can fit the trace and removes
%the points that might not be used
inNChunk = floor(length(db1Trace) / inStepLen);
db1Trace = db1Trace(1:inNChunk * inStepLen); %Removes unused points
inNChunk = inNChunk - (inNStep - 1);

%Substract the mean of the the trace in order to remove the zeros frequency
%componant
db1Trace = db1Trace - nanmean(db1Trace);

%Create a matrix, whose columns are samples of the signal whose length is
%windowLengthSec. Successive samples (columns) are overlapping and spaced
%by windowLength/step seconds. Each sample is then multiplied by a Hamming
%window.
in1ChunkBegIdx  = ((1:inNChunk) * inStepLen) - (inStepLen - 1);
in1ChunkIdx     = ((1:inWinLen) - 1)';
in2ChunkIdx     = repmat(in1ChunkBegIdx, inWinLen, 1) + repmat(in1ChunkIdx, 1, inNChunk);

%Convolutes with a hamming taper
db1Taper = hamming(size(in2ChunkIdx, 1), 'periodic');
% db1Taper = triang(size(in2ChunkIdx, 1));
if size(in2ChunkIdx, 2) ~= 1
    db2TaperTrace = db1Trace(in2ChunkIdx).*repmat(db1Taper , 1, size(in2ChunkIdx, 2));
else
    db2TaperTrace = db1Trace(in2ChunkIdx)'.*repmat(db1Taper, 1, size(in2ChunkIdx, 2));
end

% Computes the FFT and calculate the spectrum as the modulus of the FFT
myFFT=fft(db2TaperTrace);
myFFT=myFFT(1:size(db2TaperTrace,1)/2+1,:);
db2Power = (1/(inSampleRate*size(db2TaperTrace,1))).*abs(myFFT).^2;
% db2Power = (1/(inSampleRate*sum(db1Taper))).*abs(myFFT).^2;
db2Power(2:end-1,:) = 2*db2Power(2:end-1,:);

% %Multitaper version
% db2ChunkTrace = db1Trace(in2ChunkIdx);
% in1Freq = 0:2:120;
% db2Power = nan(length(in1Freq), size(db2ChunkTrace, 2));
% for iChunk = 1:size(db2ChunkTrace, 2)
%     db2Power(:, iChunk) = pmtm(db2ChunkTrace(:, iChunk), 1.5 , in1Freq, inSampleRate);
% end

%Gets the central index of each window
in1CntrIdx = ceil(median(in2ChunkIdx, 1));

%Plots the result if needed 
if blDoPlot
    dbTraceLenSec = length(db1Trace)/inSampleRate;
    inTopFreq = 120;
    figure
    imagesc(dbWinLenSec/2: dbWinLenSec/inNStep : dbTraceLenSec-(dbWinLenSec/2), 0:1/dbWinLenSec:inTopFreq, ...
        10*log10(db2Power(1:inTopFreq*dbWinLenSec+1,:)));
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    set(gca, 'YDir', 'normal')
    colorbar
end