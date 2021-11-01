import numpy as np
import phate
import scprep
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import umap
from sklearn.decomposition import PCA
from scipy import stats
import time
import os

# def PlotEmbed(sTROUGH, in1EmbedLabel, blLegend, blColorbar, chLegend, chLabel, chTitle, chMethod, blZScore, inN_Component, dbAlpha, sOPTION): 
def PlotEmbed(sTROUGH, in1EmbedLabel, chLabelTag, chEmbedMethod, chDataFormat, blZScore, inN_Component, chFormatImg, blRotate3D, blAddLegend, inFontSize, blDiscrete, dbAlpha, sOPTION):

    start_PlotEmbed = time.time()
    
    # Check if the folder exists
    if not os.path.exists(os.path.join(sOPTION.chOutPath,'plots')):
        os.makedirs(os.path.join(sOPTION.chOutPath,'plots'))
        
    '''
    Convert data to the wanted representation if wanted (the complex
    trough can be represented in complex or polar coordinates)
    '''
    if blZScore: 
        db2Data =stats.zscore(sTROUGH.db2Trough, axis=0, ddof=1)
    else:
        db2Data = sTROUGH.db2Trough
    
    time_stamps     = sTROUGH.in1Index
    labels          = in1EmbedLabel # Give an option for just frame 
    if labels is None:
        labels = np.arange(db2data.shape[0]) # Just frames


    #Perform embedding with chosen method
    if chEmbedMethod=='umap':
        print('Running UMAP ({}D)'.format(inN_Component))
        reducer = umap.UMAP(n_components=inN_Component, n_neighbors=5,random_state=0)
        emb_data = reducer.fit_transform(db2Data)

    elif chEmbedMethod=='pca':
        print('Running PCA ({}D)'.format(inN_Component))
        pca = PCA(n_components=inN_Component)
        emb_data = pca.fit_transform(db2Data)

    else:
        print('Running PHATE ({}D)'.format(inN_Component))
        phate_operator = phate.PHATE(n_components=inN_Component,gamma=0, random_state=0, n_jobs=-2)
        emb_data = phate_operator.fit_transform(db2Data)
        

    if '20Clusters' in chLabelTag:
        cmap = 'tab20'
    else:
        cmap = 'viridis'
    
    # Plot the embedding
    if inN_Component == 3:    
        if blRotate3D == True:
            scprep.plot.rotate_scatter3d(emb_data, c=labels, filename = sOPTION.chOutPath + '/plots/plot_embedding_' + str(inN_Component) + 'D_' + chEmbedMethod + '_' + chLabelTag + '.gif', 
                                         figsize=(10, 8), discrete = chDiscrete, fontsize = inFontSize, legend = blAddLegend, colorbar = blAddLegend, cmap=cmap, s=10, alpha=dbAlpha)
        else:
            if chFormatImg == 'eps':
                scprep.plot.scatter3d( emb_data, c=labels, filename = sOPTION.chOutPath + '/plots/plot_embedding_' + str(inN_Component) + 'D_' + chEmbedMethod + '_' + chLabelTag  + '.' + chFormatImg, 
                    figsize=(10, 8), discrete = blDiscrete, fontsize = inFontSize, legend = blAddLegend, colorbar = blLegend, cmap=cmap, s=10)
            else:
                scprep.plot.scatter3d( emb_data, c=labels, filename = sOPTION.chOutPath + '/plots/plot_embedding_' + str(inN_Component) + 'D_' + chEmbedMethod + '_' + chLabelTag  + '.' + chFormatImg,
                    dpi = 1000, discrete = blDiscrete, colorbar = blAddLegend, fontsize = inFontSize, figsize=(10, 8), cmap=cmap, s=10, alpha=dbAlpha)
    else:
        if chFormatImg == 'eps':
            scprep.plot.scatter2d(emb_data, c=labels, filename = sOPTION.chOutPath + '/plots/plot_embedding_' + str(inN_Component) + 'D_' + chEmbedMethod + '_' + chLabelTag + '.' + chFormatImg, 
                                         figsize=(10, 8), discrete = blDiscrete, fontsize = inFontSize, legend = blAddLegend, colorbar = blLegend, cmap=cmap, s=10)
        else:
            scprep.plot.scatter2d(emb_data, c=labels, filename = sOPTION.chOutPath + '/plots/plot_embedding_' + str(inN_Component) + 'D_' + chEmbedMethod + '_' + chLabelTag + '.' + chFormatImg, 
                                         dpi = 1000, discrete = blDiscrete, fontsize = inFontSize, legend = blAddLegend, colorbar = blLegend, figsize=(10, 8), cmap=cmap, s=10, alpha=dbAlpha)

    plt.close('all')
    if sOPTION.blVerbose: print("---- PLotting time: {}s seconds ---\n".format(time.time() - start_PlotEmbed))
    print('Done!')

    
    return emb_data