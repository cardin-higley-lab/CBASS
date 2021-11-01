from scipy import stats
import numpy as np

from Pipeline.Utilities.CBASS_U_EnrichmentScore import EnrichmentScore
from Pipeline.Utilities.CBASS_U_ScoreThreshold import ScoreThreshold


def PartitionTrough(sTROUGH, bl1Epoch, bl1Baseline=None, blZScore=True, inNClu=20, dbSigThrs=10^-4, inNIter=1000, sOPTION=None):
    '''
    L2 of the bout pipeline: Estimate the probability for a set of events to
    occur during a state indexed by bl1Epoch. Then find a threshold that
    best separates events having a high and a low probability of occurence.
    Events are meant to be the output of CBASS_L1_GetTrough and correspond to
    the troughs of oscilatory activity at the band of interest in a
    reference channel.

    Input -------------------------------------------------------------------

    sTROUGH:      the output of CBASS_L1_GetTrough (i.e.) a structure
                  requiring the following fields:
                  -.db1FilterBand an (1 x 2) array describing the frequency
                  band of interest i.e. [30 80] for 30 to 80 Hz.
                  -.db2Trough  a (2 * channel x trough) matrix containing the
                  hilbert transform of each channel of sREC.db2LFP filtered
                  in the band defined in sTROUGH.db1FilterBand at the trough
                  of the filtered signal in a reference channel (see
                  CBASS_L1_GetTrough for more detail)
                  -.in1Index the indices of the trough in sREC.db2LFP
    bl1Epopch:    a logical vector, containing as many elements as time
                  samples in db2LFP, indexing the state in which enriched
                  band specific activity is observed.
    bl1Baseline:  (optional) a logical vector, containing as many elements as
                  time samples in db2LFP, indexing the state in which band
                  specific activity is not observed. There should be no
                  overlap between bl1Baseline and bl1Epoch. Default is
                  ~bl1Epoch.
    blZScore:     (optional) logical specifying if trought data is to be
                  zscored for k-means partitioning. Default is true.
    inNClu:       (optional) number of cluster used for k-means partitioning.
                  Default is 20
    dbSigThrs:    (optional) threshold for the significance of the enrichment
                  in trough partition.  P-Values are computed with a binomial 
                  test. Default is 10.^-4.
    inNMaxIter:   (optional) maximum iteration used by the k-means algorithm.
                  Default is 1000.
    blVerbose:    (optional) logical setting whether processing updates
                  should be displayed in the command window. Default is true.

    Output ------------------------------------------------------------------

    sPART     	a structure containing the following fields:
                  -.db1Score a real vector giving an estimate of the
                  probability to occur in the state of interest for  each
                  trough in STROUGH based on its position in the feature
                  space
                  -.dbThreshold a scalar giving the threhsold value of
                  db1Score above witch troughs are considered part of the
                  partition
                  -.bl1Partition a boolean indexing troughs in the partition
    '''


    # Sets optional parameters if not provided
    bl1Epoch = bl1Epoch[0].flatten()
    if np.any(bl1Baseline==None) or len(bl1Baseline)==0: 
        print('bl1Baseline is being set to default')
        bl1Baseline = ~bl1Epoch
    else:
        bl1Baseline = bl1Baseline.flatten()

    # Checks that bl1Epoch and bl1Baseline are non empty and non overlapping
    if ~np.any(bl1Epoch): print('bl1Epoch does not index any time sample')
    if ~np.any(bl1Baseline): print('bl1Baseline does not index any time sample')
    if np.any(bl1Epoch & bl1Baseline): print('bl1Epoch and bl1Baseline are overlapping')

    # Formats the data
    if blZScore: 
        db2Data = stats.zscore(sTROUGH.db2Trough)
    else: 
        db2Data = sTROUGH.db2Trough

    # Computes the indices of the troughs recorded during the enriched epoch and during baseline
    bl1T_Epoch      = bl1Epoch[sTROUGH.in1Index.flatten()]
    bl1T_Baseline   = bl1Baseline[sTROUGH.in1Index.flatten()]

    # Computes the estimate of the probability of the state defined by bl1Epoch based on the shape of troughs
    db1Score = EnrichmentScore(db2Data, bl1T_Epoch, bl1T_Baseline, inNClu, dbSigThrs, inNIter, sOPTION)
    
    # Find the threshold that best separates events having low and high probabilities of occuring during bl1Epoch
    dbThreshold     = ScoreThreshold(sTROUGH.db2Trough, db1Score, sOPTION)

    # Computes a boolean of trought above threshold
    bl1Partition    = db1Score > dbThreshold

    # Store the output
    class sPART:
        pass
    sPART.db1Score      = db1Score
    sPART.dbThreshold   = dbThreshold
    sPART.bl1Partition  = bl1Partition
    
    return sPART