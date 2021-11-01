import numpy as np

def MakeEpochVector(in1OnPts, in1OffPts, inVecLen, sOPTION):
    '''
    [BL1EPOCHVECTOR] = NS_MakeEpochVector(IN1ONPTS, IN1OFFPTS, INVECLEN)
    creates a boolean vector BL1EPOCHVECTOR defining epoch out of a set of
    'on' and 'off' indices provided in the input vectors IN1ONPTS and
    IN1OFFPTS. The length of the output vector is set by the parameter
    INVECLEN.
    '''
    
    verbose = sOPTION.blVerbose
    # Ensures that the inputs are vectors
#     if ~isvector(in1OnPts) || ~isvector(in1OffPts)
#         error('in1OnPts and in1OffPts must be vectors')
#     end

#     # Ensures that in1OnPts and in1OffPts are horizontal vectors
#     if isrow(in1OnPts),     in1OnPts = in1OnPts(:)'; end
#     if isrow(in1OffPts),    in1OffPts = in1OffPts(:)'; end

    # Generates an error if the 'on' and 'off' vectors do not have the same lenght
    if len(in1OnPts) != len(in1OffPts):
        print('Mismatch between the number of ''on'' and ''off'' points')

    # Checks that there is no mismatch between nan points in the vector
    if any(np.isnan(in1OnPts) != np.isnan(in1OffPts)): print('Mismatch in NaN in ''on'' and ''off'' point vectors')

    # Removes nan
    in1OnPts[np.isnan(in1OnPts)]   = []
    in1OffPts[np.isnan(in1OffPts)]  = []
    # if verbose: print('in1OnPts: ',in1OnPts)
    # if verbose: print('in1OffPts: ',in1OffPts)

    # Checks that all off points are superior to off points
    if np.any(in1OnPts > in1OffPts): print(' ''off'' points must not preceed ''on'' points')

    # Creates of vector of boolean set to 'false' and sets epoch between matching 'on' and 'off' points to 'true'
    bl1EpochVector = np.zeros((1, inVecLen), dtype=bool)[0]
    for ii in range(len(in1OnPts)):
        # print('in1OnPts[ii], in1OffPts[ii]: ', in1OnPts[ii],in1OffPts[ii])
        bl1EpochVector[in1OnPts[ii]:in1OffPts[ii]] = True

    return bl1EpochVector