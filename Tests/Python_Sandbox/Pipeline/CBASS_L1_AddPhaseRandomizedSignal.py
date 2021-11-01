from sklearn.decomposition import PCA
import numpy as np
import time
from Pipeline.Utilities.CBASS_U_PhaseRandomize1D import PhaseRandomize1D

def AddPhaseRandomizedSignal(sREC, sOPTION):
    '''
     Part of L1 of the bout pipeline. Creates surrogate LFP having the same
     covariance matrix and the same spectral density as sREC.db2LFP.
    
     Input -------------------------------------------------------------------
    
     sREC:         a structure requiring the following fields:
                   -.db2LFP a (channel x time sample) matrix containing the
                   signal (i.e. time series) of interest. 
                   -.inSampleRate a positive number representing the sample
                   rate of the time series.
     inSeed        (optional) a non negative integer used to seed the random
                   number generator.
     
     Output ------------------------------------------------------------------
    
     sREC          the input structure with the additonal subfield:
                   -.db2LFP_Rnd a (channel x time sample) matrix of the same
                   size as sREC.db2LFP containing a phase randomized signal
                   having the same spectral density and the same covariance
                   matrix as sREC.db2LFP.
    '''

    verbose = sOPTION.blVerbose
    start_time = time.time()

    # Decomposes the LFP into orthogonal principal components
    pca = PCA(random_state = 1949,n_components=sREC.db2LFP.shape[0])
    db2LFP_Rnd = pca.fit_transform(sREC.db2LFP.T)
    if verbose: print("-- PCA: {}s seconds ---".format(time.time() - start_time))
    # if verbose: print('db2LFP_Rnd.shape: ',db2LFP_Rnd.shape)
    # if verbose: print('db2LFP_Rnd[:5,:5]: ',db2LFP_Rnd[:5,:5])
    
    # Randomizes the phase of each principal component
    start_phase = time.time()
    db2LFP_Rnd = PhaseRandomize1D(db2LFP_Rnd, sOPTION, inDim = 0)
    if verbose: print("-- Randomize1D: {}s seconds ---".format(time.time() - start_phase))
    # if verbose: print('(randomized) db2LFP_Rnd[:5,:5]: ',db2LFP_Rnd[:5,:5])
    
    # Recomposes the signal and stores it in the output structure
    sREC.db2LFP_Rnd = pca.transform(db2LFP_Rnd).T

    return sREC