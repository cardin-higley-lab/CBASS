import numpy as np
import phate
import scprep
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import umap
from sklearn.decomposition import PCA
from scipy import stats

import os

def PlotEmbed(sTROUGH, in1EmbedLabel, blLegend, blColorbar, chLegend, chLabel, chTitle, chMethod, blZScore, inN_Component, dbAlpha, sOPTION): 

    # Export variables to temp for test
    verbose = sOPTION.blVerbose
    
    # Check if the folder exists
    if not os.path.exists('./plots'):
        os.makedirs('./plots')
        
    '''
    Convert data to the wanted representation if wanted (the complex
    trough can be represented in complex or polar coordinates)
    '''
    if blZScore: 
        db2Data =stats.zscore(sTROUGH.db2Trough, axis=0, ddof=1)
    else:
        db2Data = sTROUGH.db2Trough
        
    # inNChan = int(db2Data.shape[1]/2)
    # db2Data    = np.concatenate((np.abs(db2Data[:, :inNChan] + 1j * db2Data[:, inNChan:]), 
    #                              np.angle(db2Data[:, :inNChan] + 1j * db2Data[:, inNChan:])),axis=1)
    
    # if verbose: 
    #     print('db2Data.shape: ',db2Data.shape)
    #     print('db2Data[:10,:10]: ',db2Data[:10,:10])
    
    time_stamps     = sTROUGH.in1Index
    labels          = in1EmbedLabel # Give an option for just frame 
    if labels is None:
        labels = np.arange(db2data.shape[0]) # Just frames


    #Perform embedding with chosen method
    if chMethod=='umap':
        print('Running UMAP ({}D)'.format(inN_Component))
        reducer = umap.UMAP(n_components=inN_Component, n_neighbors=5,random_state=0)
        emb_data = reducer.fit_transform(db2Data)

    elif chMethod=='pca':
        print('Running PCA ({}D)'.format(inN_Component))
        pca = PCA(n_components=inN_Component)
        emb_data = pca.fit_transform(db2Data)

    else:
        print('Running PHATE ({}D)'.format(inN_Component))
        phate_operator = phate.PHATE(n_components=inN_Component,gamma=0, random_state=0, n_jobs=-2)
        emb_data = phate_operator.fit_transform(db2Data)
        
    # Plot the embedding
    if inN_Component == 3:    
        scprep.plot.scatter3d(emb_data, c=labels, filename = './plots/plot_embedding_' + str(inN_Component) + 'D_' + chMethod + '_' + chLabel + '.png',
            dpi = 500, colorbar = blLegend, fontsize = 12, figsize=(10, 8), s=10, alpha=dbAlpha)
    else:
        scprep.plot.scatter2d(emb_data, c=labels, filename = './plots/plot_embedding_' + str(inN_Component) + 'D_' + chMethod + '_' + chLabel + '.png',
                              dpi = 500, fontsize = 12, legend = blLegend, colorbar = blColorbar, figsize=(10, 8), s=10, alpha=dbAlpha, 
                              xlabel='C0', ylabel='C1', legend_title=chLegend, title=chTitle)

    plt.close('all')
    
    return emb_data