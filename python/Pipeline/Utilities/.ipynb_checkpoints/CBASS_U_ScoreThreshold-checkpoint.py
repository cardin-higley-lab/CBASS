import numpy as np
from scipy.optimize import fmin

def ScoreThreshold(db2Data, db1Score, sOPTION):
    '''
    CBASS utility. Calculates a threshold for the enrichment score DB1SCORE
    returned by the function CBASS_U_EnrichmentScore. The function returns 
    the threshold that maximizes the normalized mahalanobis distances between 
    the event above and under the threshold. The normalized distance can be 
    thought of as multidimensional analog of the t-statistics. The threshold 
    chosen optimizes the distance while taking sampling variability into 
    account. This version of the function uses fminsearch.

    Input -------------------------------------------------------------------

    db2Data:      a matrix (observation x parameter) representing the a set
                  of observations in a parameter space
    db1Score:     a score valued between 0 and 1 assigned to each observation
    '''

    # Performs a search
    def MDNORM(dbThr):
        
        # Sets the distance to infinity if it is out of bound
        if dbThr <= 0 or dbThr >= 1: 
            db_MDNorm = -np.inf
        else:
            # Gets troughs above threshold and their rate
            bl1TSelIdx  = db1Score > dbThr

            # Computes the coordinates of the centroids of trough above and under enrichment score threshold
            db1CntrSup  = np.mean(db2Data[bl1TSelIdx.flatten(),:])
            db1CntrMin  = np.mean(db2Data[~bl1TSelIdx.flatten(),:])

            # Computes the euclidean and mahalanobis distance between the centroids
            db_MahalD   = MahalanobisD(db2Data, db1CntrSup - db1CntrMin)
            db_MDNorm   = -db_MahalD / np.sqrt(1/np.sum(bl1TSelIdx) + 1/np.sum(~bl1TSelIdx))
            
        return db_MDNorm

    def MahalanobisD(db2RefPop, db2Obs):
        # DB1MAHAD = MahalanobisD(DB2REFPOP, DB2OBS)
        # Returns DB1MAHAD: the Mahalanobis distances of a set of observation
        # DB2OBS relative to a set of reference observation DB2REFPOP.

        db2CCov     = np.ma.cov(db2RefPop.T) # Covariance matrix of the cluster
        db1CMu      = np.nanmean(db2RefPop, axis=0) #Centroid of the cluster

        db2ObsCnt   = db2Obs - db1CMu #Centered coordintates
        db1MahaD    = np.sqrt(np.nansum(np.dot(db2ObsCnt,np.linalg.inv(db2CCov))*db2ObsCnt)) #Squared Mahalonobis distances

        return db1MahaD

    dbThreshold = fmin(MDNORM, np.median(db1Score), maxiter=50, disp=sOPTION.blVerbose)
    return dbThreshold
    # - Utility function ------------------------------------------------------

