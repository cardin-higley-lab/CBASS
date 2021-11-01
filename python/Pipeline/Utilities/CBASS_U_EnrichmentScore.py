import warnings
import numpy as np
# from scipy import stats
# from random import sample
from scipy.spatial.distance import cdist
from scipy.stats import norm, zscore

def EnrichmentScore(db2Data, bl1State, bl1Baseline, inNClu, dbSigThrs, inNItr, sOPTION):
    '''
    CBASS utility. Calculates an enrichment score for each observations (i.e.
    rows) in the data matrix DB2DATA. The score represent how likely an
    observation is to fall into a region that has more observation labelled
    by the boolean BL1STATE than by the boolean BL1BASELINE when the data
    manifold is partitionned into INNCLU regions. The score is estimated by
    performing INNITR random partition of the data, and testing for a higher
    rate of occurence of BL1STATE within each region using a binomial
    test. Partitions are performed using a method analogous to the first
    step of the k-means algorithm. INNCLU centers are drawn from the
    observation at random and observations are assigned to their closest
    center.

    Input -------------------------------------------------------------------

    db2Data:      a matrix (observation x parameter) representing the a set
               of observations in a parameter space
    bl1State:     a boolean vector indexing each observation
    bl1Baseline:  (optional) a boolean vector indexing each observation and
               having no overlap with bl1State. Default is ~bl1State.
    inNClu:       (optional) the number of regions used to conmpute the
               score. Default is 20. Higher values will yield steeper
               score distributions.
    dbSigThrs:    (optional) threshold for the significance of the enrichment
               in each region.  P-Values are computed with a binomial 
               test. Default is 10.^-4.
    bl1Zscore:    (optional) logical specifying if data is to be zscored for 
               partitioning. Default is true.
    inNItr:       (optional) integer specifying how many iteration of
               the algorithm are performed. Default is 1000;
    '''
    
    # Checks argurments
    
    inNObs = db2Data.shape[0]
    if len(bl1State) != inNObs:
        print('bl1State must be a boolean vector having as many elements as there is rows in db2Data')

    bl1State = bl1State.flatten() # To ensure it is a row vector
    if len(bl1Baseline) != inNObs:
        print("bl1Baseline must be a boolean vector the same size aas bl1State. Set to default")
        bl1Baseline = bl1State
    elif np.any(bl1Baseline & bl1State):
        print('Bl1Baseline and bl1State cannot have common elements. Set to default')
        bl1Baseline = np.invert(bl1State)

    bl1Baseline = bl1Baseline.flatten()

    # Calculates the global rate
    dbRate_All = np.sum(bl1State)/np.sum(bl1State | bl1Baseline)

    # Zscores data if needed
    if sOPTION.blZScore: db2Data = zscore(db2Data, axis=0, ddof=1)

    # Initializes the count vector
    in1SigCount = np.zeros([inNObs, 1]);

#      For a number of time set by inNItr:
#      1. Draw inNClu centers from the distribution
#      2. Assign each points to its closest center to define a region
#      3. Checks whether each region is enriched
#      4. Adds a count for each point situated in an enriched region
    for iItr in range(inNItr):
        # Select intialization centroids
        db2Cntr = DrawCentroids(db2Data, inNClu);

        # Calculates the distance of each point to the centroid and assigns each point to its cluster
        db2D = cdist(db2Data, db2Cntr)
        in1Clu = np.argmin(db2D, axis=1)

        # Calculates the enrichment of each regions and its significance
        for iClu in range(inNClu):
            bl1Clu = in1Clu == iClu #bl1Clu = np.where(in1Clu == iClu)[0]
            inNObs = np.sum(bl1Clu) #inNObs = len(bl1Clu)

            dbRate = np.sum(bl1State[bl1Clu])/np.sum(bl1State[bl1Clu] | bl1Baseline[bl1Clu])
            dbXPVal = BinomialTest(inNObs, dbRate, dbRate_All)
                
            if (dbRate > dbRate_All) &  (dbXPVal < dbSigThrs):
                in1SigCount[bl1Clu] = in1SigCount[bl1Clu] + 1

    # Normalizes the count by the number of interation to obtain an enrichment score for each point
    in1Score = in1SigCount / inNItr

    return in1Score

# Utilities ---------------------------------------------------------------

# Function 1 - Pick a set of centroids spanning the volume of the space
def DrawCentroids(db2Data, inNClu):

    # Select inNCLU random observation to serve as centers
    in1Idx  = np.random.choice(np.arange(db2Data.shape[0]), inNClu, replace=False)
    db2Cntr = db2Data[in1Idx, :]

    return db2Cntr

# Function 2 - Performs a Binomial test for one sample and a test rate
def BinomialTest(inXNObs, dbXRate, dbXRate_Tot):
    dbXPVal = 2 * (1 - norm.cdf(np.abs(dbXRate-dbXRate_Tot)/np.sqrt(dbXRate_Tot * (1 - dbXRate_Tot)/ inXNObs)))
    
    return dbXPVal

   