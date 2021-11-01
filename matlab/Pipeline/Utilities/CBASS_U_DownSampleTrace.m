function db1OutputTrace = CBASS_U_DownSampleTrace(db1InputTrace, inDownSamplingFactor, varargin)
%DB1OUTPUTTRACE = CBASS_U_DownSampleTrace(DB1INPUTTRACE, INDOWNSAMPLINGFACTOR,[INTYPE])
%
%Downsamples DB1INPUTTRACE by a factor INDOWNSAMPLINGFACTOR and outputs it
%in DB1OUTPUTTRACE. Tosses the extra points if INDOWNSAMPLINGFACTOR divides
%the length of DB1INPUTTRACE with a rest. The optional INTYPE argument
%specifies the type of downsampling. If set to 0 it takes the mean of each
%chunk, if the set to -1 it takes the minimum of each chunk and if set to 1
%it takes the minimum. It takes the mean of the chunk of the INTYPE is left
%unspecified or if it is set to any other value.

narginchk(2,3)
if inDownSamplingFactor == 1
    db1OutputTrace = db1InputTrace;
elseif mod(inDownSamplingFactor, 1) == 0 && inDownSamplingFactor > 1;
    inLength = floor(length(db1InputTrace)/inDownSamplingFactor)*inDownSamplingFactor;
    if nargin == 3
        switch varargin{1}
            case -1
                db1OutputTrace = min(reshape(db1InputTrace(1:inLength), inDownSamplingFactor, inLength/inDownSamplingFactor));
            case 1
                db1OutputTrace = max(reshape(db1InputTrace(1:inLength), inDownSamplingFactor, inLength/inDownSamplingFactor));
            otherwise
                db1OutputTrace = mean(reshape(db1InputTrace(1:inLength), inDownSamplingFactor, inLength/inDownSamplingFactor));
        end
    else
        db1OutputTrace = mean(reshape(db1InputTrace(1:inLength), inDownSamplingFactor, inLength/inDownSamplingFactor));
    end
else
    error('inDownSamplingFactor must be a positive integer')
end