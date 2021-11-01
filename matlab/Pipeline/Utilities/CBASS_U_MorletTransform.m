function cp2MorletSTR = CBASS_U_MorletTransform(db1Trace, inSampleRate, varargin)
%CP2MORLETSTR = CBASS_U_MorletTransform(DB1TRACE, INSAMPLERATE, [DB1FREQ])
%returns CP2MORLETSTR: the spectro temporal representation of DB1TRACE,
%computed with a Morlet Wavelet. The core frequency of the wavelet is set
%by DB1FREQ which must be a one dimensional vector of frequencies (in Hz).
%The default is 2:2:120; The sample rate for DB1TRACE must be specified in
%INSAMPLERATE.

%Checks the number of arguments
narginchk(2, 3)

%Sets the frequency range for analysis
if nargin == 3
    db1Freq = varargin{1};
else 
    db1Freq = 2:2:120;
end

%Makes sure that db1Trace is a row vector
db1Trace = db1Trace(:)';

%Defines the chunk size (processing by chunck spares memory)
inNumCycle  = 9;
inChunkSize = 1000000;
inOLSize    = round(inNumCycle * inSampleRate / db1Freq(1)); %OL = overlap between chunks

%Sets the parameters of the transform
fc = centfrq('cmor1-2'); 
scales = fc./(db1Freq*(1/inSampleRate));

%Does the processing in one go if their is only one chunk and in chunks
%otherwise
inNChunk    = ceil(length(db1Trace)/inChunkSize); %%%%%%%%%%%%%%%%%%%%%%%%%%
if inNChunk > 1
    %Computes the beginning of each chunk
    in1ChunkBeg = (1:inChunkSize:inNChunk*inChunkSize);
    in1ChunkEnd = [(in1ChunkBeg(2:end)-1) length(db1Trace)];
    
    %Pads the trace on the edges with the a reverted version of the
    %beginning and the end of the trace and adjust indices for the
    %beginning and the end of the trace accordingly
    db1Trace2    = [db1Trace(inOLSize + 1:-1:2) db1Trace db1Trace(end - 1:-1:end -inOLSize)];
    in1PadChunkBeg = in1ChunkBeg;
    in1PadChunkEnd = in1ChunkEnd + 2*inOLSize;
    
    %Computes the actual transform chunk by chunk and removes the overlap
    cp2MorletSTR = zeros(length(scales), length(db1Trace));
    cp2MorletSTR = complex(cp2MorletSTR);
    for iChunk = 1:inNChunk
        cp2ChunkSTR = cwt(db1Trace2(in1PadChunkBeg(iChunk):in1PadChunkEnd(iChunk)),scales,'cmor1-2'); %Spectro Temporal Representation
        cp2MorletSTR(:,in1ChunkBeg(iChunk):in1ChunkEnd(iChunk)) = cp2ChunkSTR(:, inOLSize + 1: end - inOLSize);
    end
else
    %Case where there is only one chunk
    cp2MorletSTR = cwt(db1Trace,scales,'cmor1-2'); %Spectro Temporal Representation
end

% % Plots a figure if needed
% figure
% image(10*log10(abs(cp2MorletSTR)), 'CDataMapping', 'Scaled');