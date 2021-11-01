import numpy as np

def GetTStampEventIndex(db1TStamp, db1EventTStamp, sOPTION):
    '''
    IN1EVENTIDX = NS_GetTStampEventIndex(DB1TSTAMP, DB1EVENTTSTAMP)

    Finds the indices of the  continues and monotonically increasing timestamp
    vector DB1TSTAMP that are closest to the event timestamps provided in
    DB1EVENTTSTAMP. Returns NaN for any event timestamps situated outside the
    range of values of DB1TSTAMP.

    '''

    # Checks if the values of db1TStamp are monotonically increasing
    verbose = sOPTION.blVerbose
    db1TStamp = db1TStamp.flatten()
    if np.any(np.diff(db1TStamp) <= 0):
        print('The values of db1TStamps should be monotonically increasing. Rescue mode enabled.')
        blRescueMode = True
        
    else: blRescueMode = False   

    # Sorts the events by timestamps
    in1SortIdx   = np.argsort(db1EventTStamp)
    # if verbose: print('in1SortIdx: ',in1SortIdx)
    db1EventTStamp = db1EventTStamp[in1SortIdx]
    # if verbose: print('db1EventTStamp: ',db1EventTStamp)

    # Gets the interval in the time stamp trace
    # db1TSInterval = np.array([np.diff(db1TStamp)/2, np.median(np.diff(db1TStamp)/2)])[0]
    dbTSInterval = np.median(np.diff(db1TStamp)/2)
    # if verbose: print('dbTSInterval: ',dbTSInterval)

    # Finds the closest index to each 
    in1EventIdx = np.empty_like(db1EventTStamp, dtype=np.int32)
    # if verbose: print('in1EventIdx.shape: ',in1EventIdx.shape)
    TSIdx = 0
    TStart = np.sum(db1EventTStamp < db1TStamp[0] - dbTSInterval) # Starts with the first index that is within the trace
    # if verbose: print('TStart: ',TStart)

    for jj in range(TStart,len(db1EventTStamp)):
        blMatch = False
        while blMatch == False and TSIdx <= len(db1TStamp):
            # if db1EventTStamp[jj] <= db1TStamp[TSIdx] + db1TSInterval:
            if np.abs(db1TStamp[TSIdx] - db1EventTStamp[jj]) <= dbTSInterval:
                blMatch = True
                in1EventIdx[jj] = TSIdx
                # if verbose: print('in1EventIdx[{}]: {}'.format(jj,in1EventIdx[jj]))
            else:
                TSIdx = TSIdx + 1
    
    return in1EventIdx