# Import general python packages
import time
import h5py, os
import numpy as np


# Import CBASS functions
from Pipeline.Utilities.PCP_U_NormalizeLFP import PCP_U_NormalizeLFP
from Pipeline.Utilities.CBASS_U_GetTStampEventIndex import GetTStampEventIndex
from Pipeline.Utilities.CBASS_U_MakeEpochVector import MakeEpochVector

# def TestFunc():

def LoadData(chDir, sOPTION):
    '''
    L0 of the bout pipeline. Loads the data and return them into a structure
    sREC.
    '''

    if sOPTION.blVerbose: print('---->> Loading data ... ')
    start_Loading = time.time()
    # Check the input
    # Check that chDir is a directory
    if not os.path.exists(chDir):
        print('{} is not a directory'.format(chDir))
        
    #Check that the file for visual stimulation meta structure exists
    chDPSFile = 'PCP_L1_DetectPresentationSet.mat'
    if not os.path.exists(os.path.join(chDir, chDPSFile)): 
        print('{} not found'.format(chDPSFile))
        blPres = False
    else: blPres = True

    # Check that the file for ephys meta structure exists
    chMDSFile = 'PCP_L4_MakeMetaDataStructure.mat';
    if not os.path.exists(os.path.join(chDir, chMDSFile)):
        print('{} not found'.format(chMDSFile))

    # Load the input
    if blPres: 
        sINPUT_1 = h5py.File(os.path.join(chDir, chDPSFile), 'r') 
        sINPUT_1 = sINPUT_1['sCFG']
    sINPUT_2 = h5py.File(os.path.join(chDir, chMDSFile), 'r') 
    sINPUT_2 = sINPUT_2['sCFG']

    # Get the overal sample rate
    inSampleRate = sINPUT_2['sL4MMDS']['inWorkSampleRate'][0][0];

    # Extract the LFP
#     db2LFP = sINPUT_2.sL4MMDS.db2LFP; # Exctracts the LFP matrix
#     db2LFP = PCP_U_NormalizeLFP(db2LFP, inSampleRate); # Filters and z-scores the LFP 
#     db2LFP = db2LFP(2:end, :); # Removes the first channel of the LFP which is set as the reference. We are left with 15 channels.
    
    db2LFP = sINPUT_2['sL4MMDS']['db2LFP'][()].T # Exctracts the LFP matrix
    db2LFP = PCP_U_NormalizeLFP(db2LFP, inSampleRate, sOPTION) # Filters and z-scores the LFP
    db2LFP = db2LFP[1:, :] # Removes the first channel of the LFP which is set as the reference. We are left with 15 channels.

    # Get the indices of when a visual presentation was shown
    db1TStamps_I2 = sINPUT_2['sL4MMDS']['db1TStamps'][()].flatten()
    db1PresOnTStamp_I1 = sINPUT_1['sL1DP']['db1PresOnTStamp'][()].flatten()
    db1PresOffTStamp_I1 = sINPUT_1['sL1DP']['db1PresOffTStamp'][()].flatten()
    # print('blPres: ',blPres)
    # print('len(db1TStamps_I2): ',len(db1TStamps_I2))
    if blPres:
        in1PresOnIdx    = GetTStampEventIndex(db1TStamps_I2, db1PresOnTStamp_I1, sOPTION)
        in1PresOffIdx   = GetTStampEventIndex(db1TStamps_I2, db1PresOffTStamp_I1,sOPTION)
        bl1Pres         = MakeEpochVector(in1PresOnIdx, in1PresOffIdx, len(db1TStamps_I2), sOPTION)
        # print('in1PresOnIdx: ',in1PresOnIdx)
        # print('in1PresOffIdx: ',in1PresOffIdx)
        # print('np.where(bl1Pres==True)[0]: ',np.where(bl1Pres==True)[0])

    else:
        bl1Pres = np.zeros(1, db2LFP.shape[1]);

    # Get the indices of when the mouse is running or whisking
    bl1Run          = sINPUT_2['sL4MMDS']['bl1WheelOn'][()].flatten().astype(bool)

#     # Compute the indices of exclusive conditions
#     bl1RunOnly 		= bl1Run & ~bl1Stim & ~bl1LickBout & ~bl1LickExcl
#     bl1StimOnly 	= ~bl1Run & bl1Stim & ~bl1LickBout & ~bl1LickExcl
#     bl1LickOnly   = ~bl1Run & ~bl1Stim & bl1LickBout & ~bl1LickExcl
#     bl1Quiet		= ~bl1Run & ~bl1Stim & ~bl1LickBout & ~bl1LickExcl
    
    class sREC:
        pass
    sREC.inSampleRate   = inSampleRate
    sREC.db2LFP         = db2LFP
    sREC.bl1Pres        = bl1Pres
    sREC.bl1Run         = bl1Run
    
    if sOPTION.blVerbose: print("---- Total processing time: {}s seconds ---\n".format(time.time() - start_Loading))
    
    return sREC