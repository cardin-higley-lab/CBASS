import time
import numpy as np
from scipy import fftpack

# The wavelet functions
import pywt
from scipy import signal
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
    
    verbose = sOPTION.blVerbose
    
    # Creates a Label if not provided
    if chLabel==None: chLabel = str(db1FilterBand[0]) + '-' + str(db1FilterBand[1]) + 'Hz' 
    if chDataFormat==None or not chDataFormat in ['polar', 'complex'] :
        pritn('{} is not a valid method set to UMAP'.format(chDataFormat))
        chDataFormat = 'complex'


    # start = time.time()
    # Filtering with Butterworth
    [b,a] = butter(2, 2*db1FilterBand / inSampleRate, 'bandpass') #Nth-order digital or analog Butterworth filter and return the filter coefficients.
    # if verbose: 
    #     print('Window: ', 2 * db1FilterBand / inSampleRate)
    #     print('Coefficients: b:{}, \na:{}'.format(b,a))

    # Filter
    if not sOPTION.TransformMETHOD == 'wavelet': # If I am using wavelet, don't do the filtering
        start_filtering = time.time()
        print('Using ', sOPTION.FilterMETHOD)
        if sOPTION.FilterMETHOD == 'lfilter':
            db2LFP_filt = lfilter(b, a, db2LFP).T # 1D filter. This matches: filter(B, A, db2LFP');
        elif sOPTION.FilterMETHOD == 'filtfilt':
            # db2LFP_filt = filtfilt(b, a, db2LFP, method="gust").T # Apply a digital filter forward and backward to a signal.
            db2LFP_filt = 2*filtfilt(b, a, db2LFP.T, padlen=1,axis=0) # 2x to match matlab's output
        if verbose: print('--Time for fitering: {}'.format(time.time()-start_filtering))

    # if verbose:
    #     print('db2LFP_filt.shape: {}'.format(db2LFP_filt.shape))
    #     print(db2LFP_filt[:5,:5])
        

    ## Hilbert transform, amplitude and phase
    print('Using ', sOPTION.TransformMETHOD)
    start_transform = time.time()
    if sOPTION.TransformMETHOD == 'hilbert':
        # db2_Hilbert = hilbert(db2LFP_filt, axis=0) #hilbert(db2LFP_filt.T).T # Old (original), try: hilbert(db2LFP_filt, axis=0)
        FastHilbert = lambda x: hilbert(x, fftpack.next_fast_len(len(x)), axis=0)[:len(x)]
        db2_Hilbert = FastHilbert(db2LFP_filt) # Old (original)
        # print('db2_Hilbert.shape: ',db2_Hilbert.shape)
        # print('db2_Hilbert[:10,:10]: ',db2_Hilbert[:10,:10])
        # FastHilbert = lambda x: signal.hilbert(x, fftpack.next_fast_len(len(x)))[:len(x)]
        # db2_Hilbert = FastHilbert(db2LFP_filt.T).T # Old (original)
        # Had to transpose since the Hilbert is computed along the last axis (ie, (15, 5249484)). The second transpose puts it back in the right orientation ((5249484, 15))

    elif sOPTION.TransformMETHOD == 'fft2': # Not functional yet 
        # print('len(db2LFP_filt): ',len(db2LFP_filt))
        db2_Hilbert = 2*np.fft.fftshift(np.fft.fft2(db2LFP_filt))/len(db2LFP_filt)

    elif sOPTION.TransformMETHOD == 'fft': # Not functional yet 
        # print('len(db2LFP_filt): ',len(db2LFP_filt))
        db2_Hilbert = 2 * np.fft.fftshift(np.fft.fft(db2LFP_filt))/len(db2LFP_filt)

    elif sOPTION.TransformMETHOD == 'wavelet':
        wav = pywt.ContinuousWavelet('cmor1.5-1.0') # Define the wavelet
        # if verbose:
        #     print('Central freq: ',pywt.central_frequency(wav, precision=10))
        #     print('In Hz: {}Hz '.format(pywt.scale2frequency(wav, 45, precision=8) * inSampleRate))

        #Compute the transform
        db2_Hilbert = []
        for idx in range(db2LFP.shape[0]):
            arr, _=pywt.cwt(db2LFP[idx,:],45,wav, sampling_period = 1/inSampleRate)
            db2_Hilbert.append(arr)
        db2_Hilbert = np.concatenate(db2_Hilbert,axis=0).T # Gets it in the usual format
        
    # Compute the amplitude and phase of the transformed signal
    db1_Amp = np.abs(db2_Hilbert)
    db1_Phase = np.angle(db2_Hilbert)
    if verbose: print('--Time for transformation: {}'.format(time.time()-start_transform))
    
    # if verbose:
    #     print('db2LFP_filt.T.shape: ',db2LFP_filt.T.shape)
    #     print('db2_Hilbert.shape: ',db2_Hilbert.shape)
    #     # print('db2_Hilbert[:10,:10]: ',db2_Hilbert[:10,:10])
    #     print('db1_Amp[:10,:10]: ',db1_Amp[:10,:10])
    #     print('db1_Phase[:10,:10]: ',db1_Phase[:10,:10])
    #     print('db1_Amp.shape: ',db1_Amp.shape)
    #     print('db1_Phase.shape: ',db1_Phase.shape)
        
        
    ## Finds the indices of troughs on the reference channels
    in1Index = np.where((db1_Phase[:-1, sOPTION.inRefChan] > 0) & (db1_Phase[1:, sOPTION.inRefChan] < 0))[0]
    # if verbose:
    #     print('in1Index: ',in1Index)
    #     print(len(in1Index))
    
    '''
    Format the troughs so that each row is a motif and col 1:15 is the
    real part and col 16:30 the imaginary part of the hilbert transform
    of the filtered signal
    '''
    # if verbose:
    #     print('np.real(db2_Hilbert[bl1RefTrough, :]).shape: ',np.real(db2_Hilbert[in1Index, :]).shape)
    #     print('np.imag(db2_Hilbert[bl1RefTrough, :]).shape: ',np.imag(db2_Hilbert[in1Index, :]).shape)
    #     print('np.real(db2_Hilbert[bl1RefTrough, :]): ',np.real(db2_Hilbert[in1Index, :]))
    #     print('np.imag(db2_Hilbert[bl1RefTrough, :]): ',np.imag(db2_Hilbert[in1Index, :]))
        
    if sOPTION.chDataFormat=='complex':
        db2Trough = np.concatenate((np.real(db2_Hilbert[in1Index, :]), np.imag(db2_Hilbert[in1Index, :])),axis=1);
    else:
        db2Trough       = np.concatenate((db1_Amp[in1Index, :], db1_Phase[in1Index, :]), axis=1)
    # if verbose: 
    #     print('db2Trough.shape: ',db2Trough.shape)
    #     print(db2Trough[:10,:10])
        

    class sTROUGH:
        pass
    sTROUGH.chLabel         = chLabel;
    sTROUGH.db1FilterBand   = db1FilterBand;
    sTROUGH.db2Trough       = db2Trough;
    sTROUGH.in1Index        = in1Index;
    
    # print('Time to compute troughs: {}'.format(time.time()-start))
        
    return sTROUGH