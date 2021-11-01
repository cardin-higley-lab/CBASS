import time
import numpy as np
from scipy import fftpack

# The wavelet functions
import pywt
from scipy.signal import butter, lfilter, hilbert, filtfilt, sosfiltfilt

def GetTrough(db2LFP, inSampleRate, db1FilterBand, inRefChan, chLabel=None, chDataFormat=None, sOPTION=None):
    '''
     L1 of the bout pipeline: Identifies events in the activity band defined
     by db1FilterBand. Events correspond to the trougths of oscillatory
     activity at the band of interest in a reference channel inRefChan. They
     are represented as the Hilbert transform of the band filtered activity of
     each channel. The representation can either use complex or polar
     coordinates depending on the value of the optional variable chDataFormat.

     Input -------------------------------------------------------------------

     db2LFP:           a (channel x time sample) matrix containing the signal
                       of interest - originaly meant to be a recording of the
                       Local Field Potential(LFP) but it can be any
                       multichannel time series.
     inSampleRate:     a positive number describing the sample rate of the
                       signal
     db1FilterBand:    a (1 x 2) array describing the frequency band of
                       interest i.e. [30 80] for 30 to 80 Hz.
     inRefChan:        a number specifying a reference channel. Events will be
                       aligned to the trought of the band specific activity in
                       that channel for trough identification. Default is the
                       last channel of db2LFP.
     chLabel           (optional) a charactaer array describing the band of
                       interest - (i.e. 'gamma' for [30 80Hz])
     chDataFormat:     (optional) a character array specifying the format of
                       the hilbert transforms output. Can be 'complex' or
                       'polar'. Default is 'complex'.

     Output ------------------------------------------------------------------

     sTROUGH:      a structure containing the following fields:
                   -.db1FilterBand an (1 x 2) array describing the frequency
                   band of interest i.e. [30 80] for 30 to 80 Hz.
                   -.db2Trough  a (2 * channel x trough) matrix containing the
                   hilbert transform of each channel of sREC.db2LFP filtered
                   in the band defined in sTROUGH.db1FilterBand at the trough
                   of the filtered signal in the reference channel inRefChan
                   -.in1Index the indices of the trough in sREC.db2LFP

     checks for the proper number of arguments
    '''
    
    # Creates a Label if not provided
    if chLabel==None: chLabel = str(db1FilterBand[0]) + '-' + str(db1FilterBand[1]) + 'Hz' 
    if chDataFormat==None or not chDataFormat in ['polar', 'complex'] :
        pritn('{} is not a valid method set to UMAP'.format(chDataFormat))
        chDataFormat = 'complex'

    # Filtering with Butterworth
    [b,a] = butter(2, 2*db1FilterBand / inSampleRate, 'bandpass') #Nth-order digital or analog Butterworth filter and return the filter coefficients.

    # Filter
    if not sOPTION.TransformMETHOD == 'wavelet': # If using wavelet, don't do the filtering
        start_filtering = time.time()
        if sOPTION.blVerbose: print('Using ', sOPTION.FilterMETHOD)
        if sOPTION.FilterMETHOD == 'lfilter':
            db2LFP_filt = lfilter(b, a, db2LFP).T # 1D filter. This matches: filter(B, A, db2LFP');
        elif sOPTION.FilterMETHOD == 'filtfilt':
            db2LFP_filt = 2*filtfilt(b, a, db2LFP.T, padlen=1,axis=0) # 2x to match Matlab's output
        if sOPTION.blVerbose: print('--Time for fitering: {}'.format(time.time()-start_filtering))

        

    ## Hilbert transform, amplitude and phase
    if sOPTION.blVerbose: print('Using ', sOPTION.TransformMETHOD)
    start_transform = time.time()
    if sOPTION.TransformMETHOD == 'hilbert':
        FastHilbert = lambda x: hilbert(x, fftpack.next_fast_len(len(x)), axis=0)[:len(x)]
        db2_Hilbert = FastHilbert(db2LFP_filt) # Old (original)

    elif sOPTION.TransformMETHOD == 'fft2': # Not functional yet 
        db2_Hilbert = 2*np.fft.fftshift(np.fft.fft2(db2LFP_filt))/len(db2LFP_filt)

    elif sOPTION.TransformMETHOD == 'fft': # Not functional yet 
        db2_Hilbert = 2 * np.fft.fftshift(np.fft.fft(db2LFP_filt))/len(db2LFP_filt)

    elif sOPTION.TransformMETHOD == 'wavelet':
        wav = pywt.ContinuousWavelet('cmor1.5-1.0') # Define the wavelet

        #Compute the transform
        db2_Hilbert = []
        for idx in range(db2LFP.shape[0]):
            arr, _=pywt.cwt(db2LFP[idx,:],45,wav, sampling_period = 1/inSampleRate)
            db2_Hilbert.append(arr)
        db2_Hilbert = np.concatenate(db2_Hilbert,axis=0).T # Gets it in the usual format
        
    # Compute the amplitude and phase of the transformed signal
    db1_Amp = np.abs(db2_Hilbert)
    db1_Phase = np.angle(db2_Hilbert)
    if sOPTION.blVerbose: print('--Time for transformation: {}'.format(time.time()-start_transform))
                
    ## Finds the indices of troughs on the reference channels
    in1Index = np.where((db1_Phase[:-1, sOPTION.inRefChan] > 0) & (db1_Phase[1:, sOPTION.inRefChan] < 0))[0]
    
    '''
    Format the troughs so that each row is a motif and col 1:15 is the
    real part and col 16:30 the imaginary part of the hilbert transform
    of the filtered signal
    '''
        
    if sOPTION.chDataFormat=='complex':
        db2Trough = np.concatenate((np.real(db2_Hilbert[in1Index, :]), np.imag(db2_Hilbert[in1Index, :])),axis=1)
    else:
        db2Trough       = np.concatenate((db1_Amp[in1Index, :], db1_Phase[in1Index, :]), axis=1)

    class sTROUGH:
        pass
    sTROUGH.chLabel         = chLabel
    sTROUGH.db1FilterBand   = db1FilterBand
    sTROUGH.db2Trough       = db2Trough
    sTROUGH.in1Index        = in1Index
        
    return sTROUGH