import numpy as np
from scipy.signal import butter, lfilter, hilbert, filtfilt, sosfiltfilt

def PCP_U_NormalizeLFP(db2LFP, inSampleRate, sOPTION, in1ChanSel=False):
    '''
    Utility meant to normalize the LFP and make sure that it is independant
    of the inpendance of contacts
    '''

    verbose = sOPTION.blVerbose
    # Handles variable argument
    if not in1ChanSel: 
        in1ChanSel = np.arange(0,db2LFP.shape[0])

    # Selects the channels
    db2LFP = db2LFP[in1ChanSel, :]

    # Defines a high pass filter
    dbHighBound = 1;
    [B, A] = butter(2, 2*dbHighBound / inSampleRate, 'high') #Nth-order digital or analog Butterworth filter and return the filter coefficients.

    # Filters the LFP
    db2LFP = filtfilt(B, A, db2LFP, axis=0, method="gust") # Apply a digital filter forward and backward to a signal.

    # ZScores the LFP;
    dbMu    = np.nanmean(db2LFP)
    dbSig   = np.nanstd(db2LFP)
    db2LFP = (db2LFP - dbMu)/dbSig;
                      
    return db2LFP