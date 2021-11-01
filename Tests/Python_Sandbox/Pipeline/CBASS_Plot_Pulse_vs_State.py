import numpy as np
from scipy.stats import norm
from scipy import signal
from matplotlib import pyplot as plt

def Plot_Pulse_vs_State(bl1Pulse, bl1State, db1WinSec, inSampleRate, dbConvWinSec, sOPTION, ax):
    # Bout pipeline utility to estimate if the pulses in bl1Pulse are enriched
    # during the state defined by bl1State. The function plots the average
    # instantaneous frequency of pulse around the state ON transitions.

#     # Checks for optional argument
    verbose = sOPTION.blVerbose
    dbConvWinSec = 0.5
    if verbose: print('db1WinSec: ',db1WinSec)
    # if bl1Pulse is not a logical vector assumes that it is a vector of
    # indices. Converts it to a logical vector of the same size as bl1State
    if not bl1Pulse.dtype==bool: 
        in1Pulse = bl1Pulse
        bl1Pulse = np.zeros_like(bl1State,dtype=bool)
        bl1Pulse[in1Pulse] = True

    # Calculates the frequency of bouts and events during quiet and running 
    dbFreq_StateON   = inSampleRate * np.sum(bl1Pulse & bl1State)/np.sum(bl1State)
    dbFreq_StateOFF  = inSampleRate * np.sum(bl1Pulse & ~bl1State)/np.sum(~bl1State)
    if verbose: print('dbFreq_StateON: {}, dbFreq_StateOFF: {}',dbFreq_StateON, dbFreq_StateOFF)

    # Relative increase of frequency
    dbRatio = dbFreq_StateON/dbFreq_StateOFF
    if verbose: print('dbRatio: ',dbRatio)
    # Significance of the increase
    dbPVal = 2 * (1 - norm.cdf(np.abs(dbFreq_StateON - dbFreq_StateOFF)/
                              np.sqrt((dbFreq_StateON / np.sum(bl1Pulse & bl1State)) + (dbFreq_StateOFF/ np.sum(bl1Pulse & ~bl1State)) )))
    if verbose: print('dbPVal: ',dbPVal)
    

    # Print the results on the screen
    if dbPVal > 0.0001: chP = str(dbPVal)
    else: chP = '<0.0001'
        
#     chSigString = 'OFF: ' + str(dbFreq_StateOFF) + 'fHz\tON: ' + str(dbFreq_StateON) + 'fHz\tRatio: ' + str(dbRatio) + 'f\t(p = ' + chP + ')\r'
    chSigString = 'OFF: {:.2f} Hz, ON: {:.2f} Hz; Ratio: {:.2f}; (p={})'.format(dbFreq_StateOFF, dbFreq_StateON, dbRatio, chP)
    if verbose: print('chSigString: ',chSigString)

    # Convolutes the trace by a square window to compute an instantaneous frequency
    rectwin = signal.windows.boxcar(int(dbConvWinSec * inSampleRate))
    if verbose: print('rectwin: ',rectwin)
    db1BoutTrace    = np.convolve(bl1Pulse, rectwin, 'same')/dbConvWinSec
    if verbose: print('db1BoutTrace: ',db1BoutTrace)

    # Defines the relative indices of the ETA 
    in1ETA_RelIdx   = np.arange(np.round(inSampleRate * db1WinSec[0]), np.round(inSampleRate * db1WinSec[1]))
    if verbose: 
        print('in1ETA_RelIdx.shape: ',in1ETA_RelIdx.shape)
        print('in1ETA_RelIdx: ',in1ETA_RelIdx)
        
    # Sets the event set
    in1RunON            = 1 + np.where(~bl1State[:-1] & bl1State[1:] == True)[0]
    if verbose: print('in1RunON: ',in1RunON)
    bl1Rem              = (in1RunON <= -in1ETA_RelIdx[0]) | (in1RunON >= (len(bl1State) - in1ETA_RelIdx[-1]))
    if verbose: 
        print('bl1Rem.shape: ',bl1Rem.shape)
        print('bl1Rem: ',bl1Rem)
    in1RunON[bl1Rem]    = []
    if verbose: 
        print('in1RunON.shape: ',in1RunON.shape)
        print('in1RunON: ',in1RunON)

    #  Calculates the event triggered avergae
#     ExpArray = (in1ETA_RelIdx[:,np.newaxis]+in1RunON[np.newaxis,:]).astype(int)
    db1ETA_Bout     = np.mean(db1BoutTrace[(in1ETA_RelIdx[:,np.newaxis]+in1RunON[np.newaxis,:]).astype(int)], axis=1)
    if verbose: 
        print('db1ETA_Bout.shape: ',db1ETA_Bout.shape)
        print('db1ETA_Bout: ',db1ETA_Bout)
    
    db1ETSEM_Bout   = np.std(db1BoutTrace[(in1ETA_RelIdx[:,np.newaxis]+in1RunON[np.newaxis,:]).astype(int)],axis=1)/np.sqrt(len(in1RunON))
    if verbose: 
        print('db1ETSEM_Bout.shape: ',db1ETSEM_Bout.shape)
        print('db1ETSEM_Bout: ',db1ETSEM_Bout)
        
    # Defines the time
    db1Time     = in1ETA_RelIdx/inSampleRate
    if verbose: print('db1Time: ',db1Time)

    # Plot the results
#     fig,ax = plt.subplots()
    ax.fill_between(db1Time, db1ETA_Bout + db1ETSEM_Bout,  db1ETA_Bout - db1ETSEM_Bout)#, [.5 .5 .5], 'LineStyle', 'none', 'FaceAlpha', .3); hold on
    ax.plot(db1Time, db1ETA_Bout, 'k-')
    db1YL = ax.get_ylim() 
    ax.plot(np.array([0, 0]), db1YL, 'r--')
    ax.set_ylabel('Instantaneous Frequency (Hz)'); ax.set_xlabel('Time (s)');

    return chSigString