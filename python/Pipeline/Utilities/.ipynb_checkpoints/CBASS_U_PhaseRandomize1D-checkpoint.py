import numpy as np
import time
import scipy

def PhaseRandomize1D(dbXMat, sOPTION, inDim):
    '''
    DBXPHASERANDMAT = CBASS_U_PhaseRandomize1D(DBXMAT, INDIM)
    Generate a phase randomized signal having the same spectral amplitude as
    DBXMAT along dimension INDIM. 

    Permutes the dimension so that the phase is randomized along the dimension
    of choice
    '''
    
    in1PermDim = np.arange(dbXMat.ndim)
    dbXMat = np.random.permutation(dbXMat) #If x is a multi-dimensional array, it is only shuffled along its first index.

    # Caclulates the indices
    inNSmp      = dbXMat.shape[0]
    inRndIdx    = np.arange(1,1 + np.floor((inNSmp - 1)/2),dtype=np.int16)
    if inNSmp%2 == 0: 
        inLastIdx   = int(inNSmp/2 + 1)
    else: 
        inLastIdx = []
        print('empty inLastIdx')

    # Computes the FFT and the randomized phase
    start_fft = time.time()
    dbXMat = dbXMat.T
    db2FFT      = scipy.fft.fft(dbXMat, workers=-1)
    db2Amp      = np.abs(db2FFT[:, inRndIdx]);
    db2RndAng   = np.random.randn(*db2Amp.shape)* 2 * np.pi - np.pi #
    if sOPTION.blVerbose: print("-- FFT: {}s seconds ---".format(time.time() - start_fft))

    '''
    Recompose a phase randomized LFP (for the reconstituted signal to be real,
    the first element of the fft must real the rest must be conjugate
    symectric; if the number of sample is even, element NSamp/2 + 1 must be
    real)
    '''

    db2FFT_Rnd  = np.concatenate((db2FFT[:, 0].reshape(len(db2FFT[:, 0]),-1), 
                                  db2Amp  * np.exp(1j * db2RndAng), 
                                  db2FFT[:, inLastIdx].reshape(len(db2FFT[:, inLastIdx]),-1) ,
                                  np.conj(np.fliplr(db2Amp) * np.exp(1j * np.fliplr(db2RndAng)))),axis=1)

    start_ifft = time.time()
    dbXPhaseRandMat = np.real(scipy.fft.ifft(db2FFT_Rnd, workers=-1)) #Ideally, use 'real_if_close', but its tolerance would have to be really high (tol=10e12). Thus, using 'real' as a imediate fix. Long term: find out what is happening with the accuracy of the computation. 
    if sOPTION.blVerbose: print("-- IFFT: {}s seconds ---".format(time.time() - start_ifft))
    return dbXPhaseRandMat.T