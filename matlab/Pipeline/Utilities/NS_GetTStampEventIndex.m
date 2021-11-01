function [in1EventIdx] = NS_GetTStampEventIndex(db1TStamp, db1EventTStamp)
%IN1EVENTIDX = NS_GetTStampEventIndex(DB1TSTAMP, DB1EVENTTSTAMP)
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
    warning('The values of db1TStamps should be monotonically increasing. Rescue mode enabled.')
    blRescueMode = true;
else, blRescueMode = false; end
% if length(unique(diff(db1TStamp))) ~= 1
%     error('The interval of db1TStamps is not constant')
% end
 
%Sorts the events by timestamps
[db1EventTStamp, in1SortIdx]   = sort(db1EventTStamp);
[~, in1UnsortIdx]              = sort(in1SortIdx);

% Sorts the time stamps if needed
if blRescueMode
    %Sorts the events by timestamps
    [db1TStamp, in1SortIdx2]   = sort(db1TStamp);
    [~, in1UnsortIdx2]         = sort(in1SortIdx2);
end

%Gets the interval in the time stamp trace
db1TSInterval = [diff(db1TStamp)./2; median(diff(db1TStamp)./2)];

%Finds the closest index to each 
in1EventIdx = nan(1, length(db1EventTStamp));
TSIdx = 1;
TStart = 1 + sum(db1EventTStamp < db1TStamp(1) - db1TSInterval(end)); %Starts with the first index that is within the trace

for jj = TStart:length(db1EventTStamp)
    blMatch = false;
    while blMatch == false && TSIdx <= length(db1TStamp)
        if db1EventTStamp(jj) <= db1TStamp(TSIdx) + db1TSInterval(TSIdx)
            blMatch = true;
            in1EventIdx(jj) = TSIdx;
        else
            TSIdx = TSIdx + 1;
        end
    end
end

%Unsort the timestamps
if blRescueMode
     in1EventIdx(~isnan(in1EventIdx)) = in1UnsortIdx2(in1EventIdx(~isnan(in1EventIdx)));
end
in1EventIdx = in1EventIdx(in1UnsortIdx);