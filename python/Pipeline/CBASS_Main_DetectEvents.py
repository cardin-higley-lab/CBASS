from scipy import stats 
import time
import numpy as np

from Pipeline.CBASS_L1_AddPhaseRandomizedSignal import AddPhaseRandomizedSignal
from Pipeline.CBASS_L1_GetTrough import GetTrough
from Pipeline.CBASS_L2_PartitionTrough import PartitionTrough


def Main_DetectEvents(db2LFP, inSampleRate, cBAND, cSTATE, sOPTION):
    '''
    Main function for the detection of enriched band specific activity motif.
    Chance levels of detection and significance levels are estimated by
    repeating procedures on surrogate data having the same spectral density 
    the same covariance matrix between channels. The surrogate data are
    generated with the function:  CBASS_L1_AddPhaseRandomizedSignal.
    The function makes use of 2 additional subroutines:

    CBASS_L1_GetTrough:       performs trough identification of real and
                              surrogate signals using the Hilbert transform.
    CBASS_L2_PartitionTrough: identifies a groups of troughs having a high
                              probability of occuring during the state of 
                              interest baseed on spectro-temporal dynamics.


    Input -------------------------------------------------------------------
    db2LFP:       a (channel x time sample) matrix containing the signal of
                  interest - originaly meant to be a recording of the Local 
                  Field Potential(LFP) but it can be any multichannel time 
                  series.
    inSampleRate: a positive number describing the sampling rate of the
                  signal
    cBAND:        a (1 x 2) array describing the frequency band of interest
                  i.e. [30 80] for 30 to 80 Hz OR several such arrays
                  grouped in a cell array. The analysis will be performed
                  for each element.
    cSTATE:       a logical vector, containing as many elements as time
                  samples in db2LFP, indexing the state in which enriched
                  band specific activity is observed OR a cell array of such
                  vector. If so the number of element of the cell array must
                  mach the number of element in cBAND.

    Option structure --------------------------------------------------------
    The variable sOPTION can be used to pass optional arguments. It may
    contain the following fields:

    .cBAND_LABEL:     a cell array of size matching cBAND containings labels
                      for the band of interest - (i.e. 'gamma' for [30 80Hz])
    .chDataFormat:    a character array specifying the format of the hilbert
                      transforms output. Can be 'complex' or 'polar'. Default
                      is 'complex'.
    .inRefChan:       a number specifying a reference channel. Events will be
                      aligned to the trought of the band specific activity in
                      that channel for trough identification. Default is the
                      last channel of db2LFP.
    .cBASELINE        a logical vector, containing as many elements as time
                      samples in db2LFP, indexing the state in which enriched
                      band specific activity is NOT observed OR a cell array
                      of such vector. If so the number of element of the cell
                      array must mach the number of element in cBAND.
    .blZScore:        logical specifying if trought data is to be zscored for 
                      k-means partitioning. Default is true.
    .inNClu:          number of cluster used for k-means partitioning.
                      Default is 20
    .dbSigThrs:       threshold for the significance of the enrichment in 
                      trough partition.  P-Values are computed with a 
                      binomial test. Default is 10.^-4.
    .inNIter:         iteration used by the k-means algorithm. Default is 
                      1000.
    .blVerbose:       logical setting whether processing updates should be
                      displayed in the command window. Default is true.

    Output ------------------------------------------------------------------
    sFREQ_BAND:       a structure array of the same size as cBAND containing
                      the following fields:
      .db1Band        The value of cBAND for that instance of the array
      .chBandLabel    The label of the band for that instance of the array 
                      (may be specified in sOPTION.cBAND_LABEL).
      .in1TroughIdx   the indices of the trough of oscillatory activity in
                      sREC.db2LFP used to define activity templates (Troughs
                      are generated by CBASS_L1_GetTrough).
      .db1Score       a score representing an estimate of the probability of
                      trough to occur during the state of interest based on
                      spectro-temporal features
      .db1ScoreRnd    the same score calculated on trought extracted from
                      surrogate data
      .dbP_KS         the p-value of a Kolmogorov-Smirnov test of the
                      difference between the distrubutions of score values in
                      the LFP and in the surrogate data
      .dbThreshold    a threshold giving the partition of values of db1Score
                      resulting in the most significant distance between
                      troughs based on spectro-temporal features
      .dbPValThr      p-value for the significance of the partition.
                      Calculated as the proportion of db1ScoreRnd above
                      dbThreshold
      .bl1Partition   a logical vector indexing the trough used to compute 
                      the template motif
      .bl1Event       a boolean, having as many elements as time samples in
                      db2LFP, indexing significant event of state
                      enriched band specific activity. <---  FINAL OUTPUT

    Optional Output ---------------------------------------------------------
    cTROUGH           a cell array containing the output of the subroutine
                      CBASS_L1_GetTrough for each frequency band in cBAND.
                      Each instance of cTROUGH is a structure sTROUGH
                      containing the following fields:
      .db1FilterBand  an (1 x 2) array describing the frequency
                      band of interest i.e. [30 80] for 30 to 80 Hz.
      .db2Trough      a (2 * channel x trough) matrix containing the
                      hilbert transform of each channel of sREC.db2LFP filtered
                      in the band defined in sTROUGH.db1FilterBand at the trough
                      of the filtered signal in the reference channel inRefChan
      .in1Index       the indices of the trough in sREC.db2LFP

    cTRGH_RND         a cell array containing the output of the subroutine
                      CBASS_L1_GetTrough for each frequency band in cBAND
                      applied to the surrogate signal generated by
                      CBASS_L1_AddPhaseRandomizedSignal. Each instance of
                      cTRGH_RND is a structure sTRGH_RND containing the same
                      fields as the elements of cTROUGH (see above)

    '''


    # Checks for non optional arguments ----------------------------------------
    cBAND, cSTATE, inNChan = CheckArg(db2LFP, inSampleRate, cBAND, cSTATE) # Function at the end of the script

    # Deals with the option structure -----------------------------------------
    sOPTION = CheckOption(sOPTION, cBAND, cSTATE, inNChan) # Function at the end for the script

    # Computes ----------------------------------------------------------------
        # Formats the LFP and computes the phase randomized signal (this will be
        # used for chance level estimation and statistical testing)
    start_AddPhaseRandomizedSignal = time.time()
    if sOPTION.blVerbose: print('---->> Compute phase randomized signal ... ')
    class sREC:
        pass
    sREC.db2LFP         = db2LFP
    sREC.inSampleRate   = inSampleRate
    sREC                = AddPhaseRandomizedSignal(sREC, sOPTION)
    if sOPTION.blVerbose: print("---->> Total processing time: {}s seconds ---".format(time.time() - start_AddPhaseRandomizedSignal))

    if sOPTION.blVerbose: print('Done formatting the LFP and computing the phase randomized signal \n')

    # Initializes the output structure
    inNBnd = len(cBAND)
    class sFREQ_BAND_struct:
        pass
    sFREQ_BAND = [sFREQ_BAND_struct() for i in range(inNBnd+1)]
    class cTROUGH_struct:
        pass
    cTROUGH = [cTROUGH_struct() for i in range(inNBnd+1)]
    class cTRGH_RND_struct:
        pass
    cTRGH_RND = [ cTRGH_RND_struct() for i in range(inNBnd+1)]


    # Loops through bands of interest 
    for iBnd in range(inNBnd):
        # Choose the state of interest
        if sOPTION.blVerbose: print('\n------ {} ---------------------------\n'.format(sOPTION.cBAND_LABEL[iBnd]))
        bl1State    = cSTATE[iBnd]

        # Extracts trougths for the real and surrogates signals
        if sOPTION.blVerbose: print('---->> Extract hilbert troughs ... ')
        
        start_GetTrough = time.time()
        if sOPTION.blVerbose: print('\n---- Real signal ... ')
        sTROUGH     = GetTrough(sREC.db2LFP, inSampleRate, cBAND[iBnd], sOPTION.inRefChan, sOPTION.cBAND_LABEL[iBnd], sOPTION.chDataFormat, sOPTION)
        if sOPTION.blVerbose: print("---- Total processing time: {}s seconds ---\n".format(time.time() - start_GetTrough))
        
        start_GetTrough = time.time()
        if sOPTION.blVerbose: print('---- Surrogate signal ... ')
        sTRGH_RND   = GetTrough(sREC.db2LFP_Rnd, inSampleRate, cBAND[iBnd], sOPTION.inRefChan, sOPTION.cBAND_LABEL[iBnd], sOPTION.chDataFormat, sOPTION)
        if sOPTION.blVerbose: print("---- Total processing time: {}s seconds ---\n".format(time.time() - start_GetTrough))
        
     
        if sOPTION.blVerbose: print('Done with extraction \n')
        # print('bl1State.shape: ',bl1State.shape)

        # Partitions troughs
        if sOPTION.blVerbose: print('---->> Partition troughs ... ')
        start_PartitionTrough = time.time()
        if sOPTION.blVerbose: print('\n---- Real signal ... ')
        sPART       = PartitionTrough(sTROUGH, bl1State, sOPTION.cBASELINE[iBnd], sOPTION.blZScore, sOPTION.inNClu, sOPTION.dbSigThrs, sOPTION.inNIter, sOPTION)
        if sOPTION.blVerbose: print("---- Total processing time: {}s seconds ---\n".format(time.time() - start_PartitionTrough))
        
        start_PartitionTrough = time.time()
        if sOPTION.blVerbose: print('---- Surrogate signal ... ')
        sPRT_RND    = PartitionTrough(sTRGH_RND, bl1State, sOPTION.cBASELINE[iBnd], sOPTION.blZScore, sOPTION.inNClu, sOPTION.dbSigThrs, sOPTION.inNIter, sOPTION)
        if sOPTION.blVerbose: 
            print("---- Total processing time: {}s seconds ---\n".format(time.time() - start_PartitionTrough))
            print('------ Done processing {} ----------------\n'.format(sOPTION.cBAND_LABEL[iBnd]))

        # Computes control variable (prints them if wanted)
        _, dbP_KS = stats.ks_2samp(sPART.db1Score.flatten(), sPRT_RND.db1Score.flatten()) #two-sample Kolmogorov-Smirnov
        dbPValThr   = np.mean(sPRT_RND.db1Score.flatten() > sPART.dbThreshold)
        if sOPTION.blVerbose:
            if dbP_KS < 0.05: 
                chKS_Sig = '' 
            else: chKS_Sig = 'NON'
            print('Score: \t\t{} SIGNIFICANT \t(p = {}, KS Test - real vs surrogate data)'.format(chKS_Sig, dbP_KS))
            if dbPValThr < 0.05:
                chSig = '' 
            else: 
                chSig = 'NON'
            print('Partition: \t{} SIGNIFICANT \t(p = {}, Fraction of surrogate troughs above threshold)\n'.format( chSig, dbPValThr))

        # Computes a boolean indexing event (the final selection of troughs)
        bl1Event    = np.zeros((1, db2LFP.shape[1]), dtype=bool) 
        bl1Event[0, sTROUGH.in1Index[sPART.bl1Partition.flatten()]] = True

        # Aggregates pulse data
        sFREQ_BAND[iBnd].db1Band        = cBAND[iBnd]
        sFREQ_BAND[iBnd].chBandLabel    = sOPTION.cBAND_LABEL[iBnd]
        sFREQ_BAND[iBnd].in1TroughIdx   = sTROUGH.in1Index
        sFREQ_BAND[iBnd].db1Score       = sPART.db1Score
        sFREQ_BAND[iBnd].db1ScoreRnd    = sPRT_RND.db1Score
        sFREQ_BAND[iBnd].dbP_KS         = dbP_KS
        sFREQ_BAND[iBnd].dbThreshold    = sPART.dbThreshold
        sFREQ_BAND[iBnd].dbPValThr      = dbPValThr
        sFREQ_BAND[iBnd].bl1Partition   = sPART.bl1Partition
        sFREQ_BAND[iBnd].bl1Event       = bl1Event

        # Aggregate intermediary steps
        cTROUGH[iBnd]   = sTROUGH
        cTRGH_RND[iBnd] = sTRGH_RND
        
    return sFREQ_BAND, cTROUGH, cTRGH_RND


    #------------------------------------------------------------------------------------------------  
def CheckArg(db2LFP, inSampleRate, cBAND, cSTATE):
    # Utility to check non optional arguments

    # Checks that the signal is a 2D matrix
    nrows, ncols = db2LFP.shape
    if not (nrows>1 and ncols>1): print('The signal must be a (channel x time sample) matrix')
    [inNChan, inNSamp] = db2LFP.shape

    # Check sample rate
    if isinstance(inSampleRate, list) or not inSampleRate>0: print('The sample rate must be a single positive number')
    
    return cBAND, cSTATE, inNChan

# ------------------------------------------------------------------------------------------------
def CheckOption(sOPTION, cBAND, cSTATE, inNChan): #, bl1Remove):
    '''
    Utility to check the options structure. Makes use of the utility
    CheckField present at the end of the script
    '''
    
    # Checks the band label
    blLblErr = False
    if not hasattr(sOPTION, 'cBAND_LABEL') or len(sOPTION.cBAND_LABEL)<len(cBAND):
        print('sOPTION.cBAND_LABEL is not valid. Setting to default\r')
        sOPTION.cBAND_LABEL = [str(cBAND[iBnd][0]) + '-' + str(cBAND[iBnd][1]) + 'Hz' for iBnd in range(len(cBAND))] 

    # Checks L1 options for formatting hilbert troughs
    if not hasattr(sOPTION, 'chDataFormat') or sOPTION.chDataFormat not in ['polar', 'complex']:
        print('sOPTION.chDataFormat is not valid. Setting to default\r')
        sOPTION.chDataFormat = ['complex' for iBnd in range(len(cBAND))] 
        
    if (not hasattr(sOPTION, 'inRefChan'): #or not isinstance(sOPTION.inRefChan, int)): 
        print(f'sOPTION.inRefChan is set to {sOPTION.inRefChan} which is not valid. Setting to default\r')
        sOPTION.inRefChan = inNChan - 1

    # Checks L2 options for the computation of filters
        # Checks baseline option argument
    if not hasattr(sOPTION, 'cBASELINE'): #or (len(sOPTION.cBASELINE) != len(bl1Remove)):        
        print('sOPTION.cBASELINE is not valid. Setting to default\r')
        sOPTION.cBASELINE = [~cSTATE[iBnd][0] for iBnd in range(len(cBAND))] 
    
    if not hasattr(sOPTION, 'blZScore') or not isinstance(sOPTION.blZScore, bool) : sOPTION.blZScore = True
    if not hasattr(sOPTION, 'inNClu') or not isinstance(sOPTION.inNClu, int) or sOPTION.inNClu <= 0 or not isinstance(sOPTION.inNClu, int): sOPTION.inNClu = 20
    if not hasattr(sOPTION, 'dbSigThrs') or sOPTION.dbSigThrs <= 0 or sOPTION.dbSigThrs > 1: sOPTION.dbSigThrs = 10^-4
    if not hasattr(sOPTION, 'inNIter') or sOPTION.inNIter <= 0 or sOPTION.inNIter: sOPTION.inNIter = 1000
    if not hasattr(sOPTION, 'blVerbose') or not isinstance(sOPTION.blVerbose, bool) : sOPTION.blVerbose = True
    
    return sOPTION
