import numpy as np

def MakeFilter(db2LFP, inSampleRate, in1EventIdx, db1Filter_WinSec):
    #Identifies pulses of activity using the average of the LFP around a selctions of the troughs

    # Gets the size of the LFP
    inNChan, inNSample = db2LFP.shape

    # Defines the window
    in1Filt_RelIdx      = np.arange(round(inSampleRate * db1Filter_WinSec[0]) , round(inSampleRate * db1Filter_WinSec[1])+1);

    # Removes events if they are too close to the edge
    bl1Rem              = np.where((in1EventIdx <= -in1Filt_RelIdx[0]) | (in1EventIdx >= inNSample - in1Filt_RelIdx[1]))[0];
    in1EventIdx = np.delete(in1EventIdx, bl1Rem)
    inNEvt = len(in1EventIdx);

    # Creates the ETA
    db3LFP_ETA = np.empty((inNChan, len(in1Filt_RelIdx), inNEvt))
    for iEvt in range(inNEvt-1):
        db3LFP_ETA[:, :, iEvt] = db2LFP[:, np.array(in1EventIdx[iEvt] + in1Filt_RelIdx).astype(int)];

    db2Filter = np.nanmean(db3LFP_ETA, axis=2)
    
    return db2Filter