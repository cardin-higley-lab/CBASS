function [db2D_Asym, db2D_Sym] = EventOverlap(cEVENT_IDX, dbWinSec, inSampleRate)
%Synopsis :
%	[DB2D_ASYM, DB2D_SYM] = EventOverlap(cEVENT_IDX, DBWINSEC, INSAMPLERATE)
%
%Compare the overlap between sets of event occuring on a time series. Each event
%set must supplied as a vector of index in the cell array CEVENT_IDX. Events are
%considered overlapping when the fall within a window of time set by DBWINSEC
%and calculated using INSAMPLERATE i.e. the sample rate of the time series.
%Distance are returned in two matrices DB2ASYM and DB2SYM. DB2ASYM is an
%asymetrical distance matrix defined as DB2ASYM(x, y) = overlap(x, y) ./
%numel(x). DB2SYM is a symetrical distance matrix defined as DB2SYM(x, y) =
%(overlap(x, y) * 2) ./ (numel(x) + numel(y));

% Verify that that all events are vectors of positives integers
if any(~cellfun(@(x) all(mod(x, 1) == 0 & x > 0), cEVENT_IDX)) || any(~cellfun(@isvector , cEVENT_IDX))
    error('All event sets must be vector of positive integers');
end

% Makes sure that all indices are row vectors
bl1Row = cellfun(@isrow, cEVENT_IDX);
cEVENT_IDX(~bl1Row) = cellfun(@(x) x', cEVENT_IDX(~bl1Row), 'UniformOutput', false);

% Get the number of event sets
inNSet = length(cEVENT_IDX);

% Get the number of event per event set and makes sure it is stored as a
% column vector
in1NEvent = cellfun(@length, cEVENT_IDX);
if isrow(in1NEvent); in1NEvent = in1NEvent'; end

% Initializes counting matrices
db2N1i2 = zeros(inNSet);

% Calculate the window in samples
inWin   = round(inSampleRate * dbWinSec);

% Sets the chunk size (we will process things in chunk to spare memory)
try
    if isunix
        [inFail, chAns] = unix('free');
        if ~inFail
            dbMaxMem    = str2double(regexprep(chAns, '.*Mem:[ ]*(\d*).*', '$1')) * 1000; % system call to get total memory. Converts kB to bytes
            dbMaxNumel  = dbMaxMem/8; %There is 8 bytes per element in arrays type double
        else, error('Total memory could not be retrieved');
        end
    else
        sMEM        = memory;
        dbMaxNumel  = sMEM.MaxPossibleArrayBytes/8; %There is 8 bytes per element in arrays type double
    end
catch
    fprintf('Could not obtain system memory specifications. Setting chunk size to 2GB\n');
    dbMaxNumel  = 2.5*10.^8; % If everything fails sets the max number of elements to 250 million. That is 2 GB.
end
dbChunkSize = dbMaxNumel/10; % Sets the chunk size as 1/10th of the max array size;

% Calculates the overlap
for iX = 1:inNSet - 1
    for iY = iX + 1:inNSet
        % Calculates the chunk size
        inNRowMax   = floor(dbChunkSize/in1NEvent(iY));
        if inNRowMax == 0, error('cEVENT{%d} exceeds max array size', iY); end
        inNChunk    = ceil(in1NEvent(iX)/inNRowMax);
        
        % initializes aggregation variables
        inSumX  = 0; % sums events in X having a one or more matches in Y
        bl1AnyY = false(1, in1NEvent(iY)); % boolean indexing whether events in Y have had matches in X in all processed chunks
        
        % Process chunk
        for iChk = 1:inNChunk
            in1ChunkIdx = (iChk - 1) * inNRowMax + 1:min(inNRowMax * iChk, in1NEvent(iX));
            bl2Overlap  = abs(cEVENT_IDX{iX}(in1ChunkIdx)' - cEVENT_IDX{iY}) < inWin;
            inSumX  = inSumX + sum(any(bl2Overlap, 2));
            bl1AnyY = bl1AnyY | any(bl2Overlap, 1);
        end
        db2N1i2(iX, iY) = inSumX;
        db2N1i2(iY, iX) = sum(bl1AnyY);
    end
end

% Creates a symetrical matrix of the overlap
db2N1i2     = db2N1i2 + diag(in1NEvent);
db2D_Asym   = db2N1i2 ./ in1NEvent; % uses native matrix expension.
db2D_Sym    = (db2N1i2 + db2N1i2')./(in1NEvent + in1NEvent');
