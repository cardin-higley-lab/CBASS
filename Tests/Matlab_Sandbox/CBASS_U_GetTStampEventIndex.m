function [in1EventIdx] = CBASS_U_GetTStampEventIndex(db1TStamp, db1EventTStamp)
%IN1EVENTIDX = CBASS_U_GetTStampEventIndex(DB1TSTAMP, DB1EVENTTSTAMP)
%
%Finds the indices of the  continues and monotonically increasing timestamp
%vector DB1TSTAMP that are closest to the event timestamps provided in
%DB1EVENTTSTAMP. Returns NaN for any event timestamps situated outside the
%range of values of DB1TSTAMP.
%
%2016-10-12 QP: Created

%Checks if the values of db1TStamp are monotonically increasing
if ~isvector(db1TStamp)
    error('db1TStamps should be a vector\r')
end
db1TStamp = db1TStamp(:);
if any(diff(db1TStamp) <= 0)
    error('The values of db1TStamps should be monotonically increasing\r')
end
% if length(unique(diff(db1TStamp))) ~= 1
%     error('The interval of db1TStamps is not constant')
% end
 
%Sorts the events by timestamps
[db1EventTStamp, in1SortIdx]   = sort(db1EventTStamp);
[~, in1UnsortIdx]              = sort(in1SortIdx);

%Gets the interval in the time stamp trace
dbTSInterval = median(diff(db1TStamp))/2;

%Finds the closest index to each 
in1EventIdx = nan(1, length(db1EventTStamp));
TSIdx = 1;
TStart = 1 + sum(db1EventTStamp < db1TStamp(1) - dbTSInterval); %Starts with the first index that is within the trace

for jj = TStart:length(db1EventTStamp)
    blMatch = false;
    while blMatch == false && TSIdx <= length(db1TStamp)
        if abs(db1TStamp(TSIdx) - db1EventTStamp(jj)) <= dbTSInterval
            blMatch = true;
            in1EventIdx(jj) = TSIdx;
        else
            TSIdx = TSIdx + 1;
        end
    end
end

%Unsort the timestamps
in1EventIdx = in1EventIdx(in1UnsortIdx);