def MultiChannelTemplateMatching(db2Signal, db2Template, blCenter, blNormalize):
'''
Synopsis: DB1SCORE = CBASS_U_MultiChannelTemplateMatching(DB2SIGNAL, DB2TEMPLATE, [BLNORM])
Returns a score DB1SCORE indicative of how well the multi-channel signal
DB2SIGNAL matches the spatio temporal template DB2TEMPLATE
Input: -DB2SIGNAL a (channel x time sample) matrix
       -DB2TEMPLATE a (channel x time sample) matrix. The number of
       channels must be the same as in DB2SIGNAL. The number of time
       sample must be inferior to half of the duration of the signal.
       -BLNORM normalizes the score by norm of db2Signal - default is true
Output:-DBSCORE a (1 x time sample) row vector. The score represents how
       well DB2SIGNAL matches DB2TEMPLATE around each time sample and is
       the dot product S . T where T represents is linearized template
       (DB2TEMPLATE(:) and S rempresents the linearized chunk of DB2SIGNAL
       matching the size of DB2TEMPLATE centered on each time points. T is
       normalized so that all its elements sum to 1.
'''
             
    verbose = sOPTION.blVerbose

    inNChan, inNSamp = db2Signal.shape[0], db2Signal.shape[1]
    if verbose: print('inNChan: {}, inNSamp: {}'.format(inNChan, inNSamp))
    inNSmpTmp = db2Template.shape[1]
    if verbose: print('inNSmpTmp: ',inNSmpTmp)

    # Pads db2Signal with zeros for ease of computation if the number of sample of the template is uneven padds in an asymetric way
    inPadBeg = np.floor(inNSmpTmp/2).astype(int)
    inPadEnd = np.ceil(inNSmpTmp/2).astype(int)
    if verbose: 
        print('np.zeros((inNChan, inPadBeg)).shape: ',np.zeros((inNChan, inPadBeg)).shape)
        print('np.zeros((inNChan, inPadEnd)).shape: ',np.zeros((inNChan, inPadEnd)).shape)
        print('db2Signal.shape: ',db2Signal.shape)
    db2Signal = np.concatenate((np.zeros((inNChan, inPadBeg)), db2Signal, np.zeros((inNChan, inPadEnd))),axis=1)
    if verbose:
        print('inPadBeg: ',inPadBeg)
        print('inPadEnd: ',inPadEnd)
        print('db2Signal.shape: ',db2Signal.shape)

    # Normalizes and transposes the template
    if verbose: 
        print('blCenter: ',blCenter)
        print('db2Template: ',db2Template)
        print('db2Template.shape: ',db2Template.shape)
        print('np.mean(db2Template): ',np.mean(db2Template))
    if blCenter: 
        db2Template = db2Template - np.mean(db2Template)
    if verbose: print('np.linalg.norm(db2Template.flatten(): ',np.linalg.norm(db2Template.flatten()))
    db2Template = db2Template.T / np.linalg.norm(db2Template.flatten())

    # Computes the centering offset if needed
    if blCenter:
        dbXCntr = np.sum(db2Signal[:, :-inNSmpTmp]);
        if verbose: print('db2Signal[:, :-inNSmpTmp].shape: ', db2Signal[:, :-inNSmpTmp].shape)
        for iSmp in range(1,inNSmpTmp): #2:inNSmpTmp
            if verbose: 
                print('db2Signal[:, iSmp:(-inNSmpTmp + iSmp)].shape:', db2Signal[:, iSmp:(-inNSmpTmp + iSmp)].shape)
                print('(np.sum(db2Signal[:, iSmp:(-inNSmpTmp + iSmp)],axis=0)).shape: ', (np.sum(db2Signal[:, iSmp:(-inNSmpTmp + iSmp)],axis=0)).shape)
            dbXCntr = dbXCntr + np.sum(db2Signal[:, iSmp:(-inNSmpTmp + iSmp)],axis=0)
            if verbose: print('dbXCntr.shape: ',dbXCntr.shape)
        dbXCntr = dbXCntr/db2Template.size;
    else:
        dbXCntr = 0;


    # Computes the score
    db1Score = np.matmul(db2Template[0, :].reshape(1,-1), (db2Signal[:, :-inNSmpTmp] - dbXCntr))
    if verbose: print('db1Score.shape: ',db1Score.shape)
    if blNormalize: db1SS = np.sum((db2Signal[:, :-inNSmpTmp] - dbXCntr)**2, axis=0)
    for iSmp in range(1,inNSmpTmp):
        db1Score = db1Score + (np.matmul(db2Template[iSmp, :], (db2Signal[:, iSmp:-inNSmpTmp+iSmp] - dbXCntr)))
        if blNormalize: 
            db1SS = db1SS + np.sum((db2Signal[:, iSmp:-inNSmpTmp+iSmp] - dbXCntr)**2,axis=0)
    if verbose: 
        print('np.max(db1Score): ',np.max(db1Score))
        print('np.max( np.sqrt(db1SS)): ', np.max(np.sqrt(db1SS)))
        print('db1SS.shape: ', db1SS.shape)
    if blNormalize: db1Score = db1Score / np.sqrt(db1SS)
    if verbose: 
        print('db1Score.shape: ',db1Score.shape)
    return db1Score
