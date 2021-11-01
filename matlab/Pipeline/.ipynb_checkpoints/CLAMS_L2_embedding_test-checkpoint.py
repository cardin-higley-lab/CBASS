# CBASS_L2_embedding_test.py;

import phate
import numpy as np
import os, time

print('\tDone with imports')

def run_embedding(cMOTIF_zscore, n_components, time_stamps, bl1Run, outPath):
    total = time.time()
    print('\tSaving bl1Run')
    np.save(os.path.join(outPath,'bl1Run'), bl1Run)
    print("\tTime to save bl1Run: {:.2f}-s".format((time.time() - total)))
    total = time.time()

    gamma=0
    phate_operator = phate.PHATE(n_components=int(n_components),gamma=gamma, n_jobs=1)
    print("\tTime to create operator: {:.2f}-s".format((time.time() - total)))
    total = time.time()

    print('\tRunning PHATE ({}D)'.format(int(n_components)))
    emb_data = phate_operator.fit_transform(cMOTIF_zscore)
    
    print("\tTime to embed: {:.2f}-min".format((time.time() - total)/60))

    np.save(os.path.join(outPath,'emb_data'), emb_data)

    labels = bl1Pres[time_stamps]
    idx_run = np.where(bl1Run[time_stamps]==1)
    
    phate.plot.rotate_scatter3d(motif_beta_3D[idx_run[0],:], c=time_stamps[idx_run[0],0], filename = outPath + "/gamma_embedding.gif", figsize=(10, 8), cmap='viridis', s=10, alpha=0.1)
    plt.close()

    
    # # bl1Pres(cTROUGH_IDX{iBnd})
    # channel = 1 # 0-> beta; 1 -> gamma
    # del bl1Pres
    # bl1Pres = gamma_vars['bl1Run'][0]

    # print(time_labels.shape)
    # labels = bl1Pres[cTROUGH_IDX_beta-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
    # idx_run = np.where(bl1Pres[cTROUGH_IDX_beta-1]==1)
    # print(idx_run[0][:10])