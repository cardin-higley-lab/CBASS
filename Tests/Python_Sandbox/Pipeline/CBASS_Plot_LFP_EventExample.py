def Plot_LFP_EventExample(db2LFP, inAnchor, db1WinSec, inSampleRate, sOPTION, ax):
    #Plotting utility of the bout pipeline to plot LFP event triggered averages

    verbose = sOPTION.blVerbose

    # Gets the number of channels and samples of the LFP
    [inNChan, inNSample] = db2LFP.shape[0], db2LFP.shape[1]
    if verbose: print('inNChan: {}, inNSample: {}'.format(inNChan, inNSample))

    # Plots exemple of the traces around bout onset
    in1Ex_RelIdx    = np.arange(round(inSampleRate * db1WinSec[0]), round(inSampleRate * db1WinSec[1]))
    if verbose: 
        print('in1Ex_RelIdx.shape: ', in1Ex_RelIdx.shape)
        print('in1Ex_RelIdx: ', in1Ex_RelIdx)
    db1Time         = in1Ex_RelIdx / inSampleRate

    # Selection the wanted chunk of the LFP and figures out its scale
    in1Sel      = (in1Ex_RelIdx + inAnchor).astype(int)
    if verbose: print('in1Sel: ',in1Sel)
    db2LFP_Sel  = db2LFP[:, in1Sel]
    db2_RangeLFP_Sel = np.amax(db2LFP_Sel, 1) - np.amin(db2LFP_Sel, 1)
    dbScale     = 0.7 * np.max(db2_RangeLFP_Sel);
    db1YL       = dbScale * np.array([-0.5, inNChan + 1.5]);

    # Plots the LFP
    for iChan in range(inNChan):
        ax.plot(db1Time, (dbScale * (inNChan - iChan -1)) + db2LFP[iChan, in1Sel], 'k')
        if verbose: print('(dbScale * (inNChan - iChan)): ',(dbScale * (inNChan - iChan)))

    ax.set_yticks(dbScale * np.arange(0,inNChan))
    if verbose: print('dbScale * np.arange(0,inNChan+1): ',dbScale * np.arange(0,inNChan+1))
    ax.set_yticklabels(np.arange(inNChan,0,-1))
    if verbose: print('np.arange(inNChan,0,-1): ',np.arange(inNChan-1,0,-1))
    ax.set_ylabel('Channel')
    ax.set_xlabel('Time (s)')

    # Enforces the limits of the y-axis
    if verbose: print('db1YL: ',db1YL)
    ax.set_ylim(db1YL)