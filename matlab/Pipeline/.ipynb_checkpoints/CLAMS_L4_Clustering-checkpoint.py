# CBASS_L2_embedding.py;


import os
import logging
import time
import numpy as np
import pandas as pd
import phate
#import torch
#import torch.nn as nn
import matplotlib.pyplot as plt
import h5py
import umap
import scprep
from sklearn.decomposition import PCA
import scipy.io as sio
from sklearn.neighbors import KernelDensity
from scipy import stats
# from sklearn.cluster import KMeans 
from sklearn import cluster

#For plotting the video
import matplotlib
from matplotlib import animation

#from kymatio import Scattering2D
import torch.optim
import math
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--tmpPath', type=str)
parser.add_argument('--outPath', type=str)
parser.add_argument('--n_components', type=int)
parser.add_argument('--cluster_algo', type=str)
parser.add_argument('--n_clusters', type=int)
parser.add_argument('--plot_clusters', type=int)
args = parser.parse_args()

# print('Done with imports')

print('args: {}'.format(args))

mat_vars = sio.loadmat(os.path.join(args.tmpPath,'tmp_vars.mat'))

# print(mat_vars)
bl1Run = mat_vars['bl1Run'][0]
time_stamps = mat_vars['time_stamps']

labels = bl1Run[time_stamps]
running=1 #1-> Running, 0-> Not running
idx_run = np.where(bl1Run[time_stamps]==running)
idx_run = idx_run[0]

emb_data = sio.loadmat(os.path.join(args.outPath,'emb_data.mat'))['emb_data'].squeeze()
idx_run = sio.loadmat(os.path.join(args.outPath,'idx_run.mat'))['idx_run']
# cluster_labels = sio.loadmat(os.path.join(args.outPath,'cluster_labels_kmeans.mat'))['cluster_labels'][0]
# print('cluster_labels[idx_run].shape: {}'.format(cluster_labels[idx_run].shape))
# print('emb_data[idx_run,:].shape: {}'.format(emb_data[idx_run,:].shape))
# print('len(cluster_labels[idx_run]): {}'.format(len(cluster_labels[idx_run])))
# print('len(emb_data[idx_run,:]): {}'.format(len(emb_data[idx_run,:])))


if args.cluster_algo=='kmeans':
    cluster_labels = cluster.KMeans(n_clusters = args.n_clusters, max_iter = 500).fit_predict(emb_data)
    sio.savemat(os.path.join(args.outPath,'cluster_labels_kmeans.mat'), {'cluster_labels': cluster_labels})
    
if args.cluster_algo=='spectral':
    cluster_labels = cluster.SpectralClustering(n_clusters=args.n_clusters,
            eigen_solver='arpack',
            affinity="nearest_neighbors").fit_predict(emb_data)
    sio.savemat(os.path.join(args.outPath,'cluster_labels_spectral.mat'), {'cluster_labels': cluster_labels})
    


if args.plot_clusters:
    scprep.plot.rotate_scatter3d(emb_data[idx_run,:].squeeze(), c=cluster_labels[idx_run].squeeze(), 
                                 filename = args.outPath + "/cluster_" + args.cluster_algo + "_gamma_embedding.gif", figsize=(10, 8), cmap='jet', s=10, alpha=0.1)
    plt.close()

# In[3]:


# Things to try:
#1) 1D unsupervised convolutional net
#2) 


# In[4]:


# # path = '/gpfs/ysm/scratch60/dietrich/ahf38/Cardin/quentin'
# # path = '/gpfs/ysm/scratch60/dietrich/ahf38/Cardin/quentin/'
# # path = '/home/ahf38/Documents/gamma_bouts/data_random'
# path = '/gpfs/ysm/home/ahf38/Documents/gamma_bouts/data'
# # file = 'gamma_vars.mat'
# # file = 'hilbert_vars.mat'
# file = 'hilbert_vars.mat'
# fullpath = os.path.join(path,file)


# In[5]:


# # gamma_vars = scipy.io.loadmat(fullpath)
# gamma_vars = sio.loadmat(fullpath)


# # In[6]:


# # gamma_vars


# # In[7]:


# # 15 channels of LFP (Channel 15 is the deepest in the cortex and Channel 1 is the most superficial)
# db2LFP = gamma_vars['db2LFP']
# print('db2LFP.shape: {}'.format(db2LFP.shape))


# # fig, ax = plt.subplots(figsize=(10, 4), dpi= 80, facecolor='w', edgecolor='k')
# # channel = 1
# # plt.plot(db2LFP[channel,:])
# # ax.set_title('LFP ')
# # ax.set_ylabel('LFP (channel {})'.format(channel))
# # ax.set_xlabel('Frame #')


# # In[77]:


# # indices of when a visual presentation was shown (Always same constrast?)
# bl1Pres = gamma_vars['bl1Pres']
# print('bl1Pres.shape: {}'.format(bl1Pres.shape))

# # fig, ax = plt.subplots(figsize=(15, 8), dpi= 80, facecolor='w', edgecolor='k')
# # plt.plot(bl1Pres[0,:])
# # ax.set_title('visual presentation')
# # ax.set_ylabel('Stimuli')
# # ax.set_xlabel('Frame #')


# # In[9]:


# # indices of when the mouse is running or whisking
# bl1Run = gamma_vars['bl1RunOnly']
# print('bl1Run.shape: {}'.format(bl1Run.shape))

# fig, ax = plt.subplots(figsize=(10, 4), dpi= 80, facecolor='w', edgecolor='k')
# plt.plot(bl1Run[0,:])
# ax.set_title('Running only')
# ax.set_ylabel('Stimuli')
# ax.set_xlabel('Frame #')


# # In[10]:


# # indices of when the mouse is running or whisking
# bl1Run = gamma_vars['bl1Run']
# print('bl1Run.shape: {}'.format(bl1Run.shape))

# fig, ax = plt.subplots(figsize=(10, 4), dpi= 80, facecolor='w', edgecolor='k')
# plt.plot(bl1Run[0,:])
# ax.set_title('Running')
# ax.set_ylabel('Stimuli')
# ax.set_xlabel('Frame #')


# # In[76]:


# # indices of when the mouse is running or whisking
# db1_Amp = gamma_vars['db1_Amp'].transpose()
# print('db1_Amp.shape: {}'.format(db1_Amp.shape))

# # fig, ax = plt.subplots(2,1,figsize=(10, 4), dpi= 80, facecolor='w', edgecolor='k')
# # ax[0].plot(db1_Amp[0,:])
# # ax[0].set_title('Amplitude')
# # ax[0].set_ylabel('Amp')
# # ax[0].set_xlabel('Frame #')

# # for i in range(db1_Amp.shape[0]):
# #     ax[1].plot(db1_Amp[i,:1000])
# # ax[1].set_title('Amplitude [zoom]')
# # ax[1].set_ylabel('Amp')
# # ax[1].set_xlabel('Frame #')

# # fig.tight_layout(pad=1.0)


# # In[12]:


# # channel = 1 # 0-> beta; 1 -> gamma
# # print('Channel: {}'.format(channel))
# # motif_beta = gamma_vars['cMOTIF_zscore'][0][channel]
# # gamma=0
# # phate_operator = phate.PHATE(n_components=2,gamma=gamma, n_jobs=-1)
# # motif_beta_2D = phate_operator.fit_transform(motif_beta)


# # In[13]:


# # bl1Pres(cTROUGH_IDX{iBnd})
# channel = 1 # 0-> beta; 1 -> gamma
# del bl1Pres
# bl1Pres = gamma_vars['bl1Run'][0]
# cTROUGH_IDX_beta = gamma_vars['cTROUGH_IDX'][0][channel]
# print(cTROUGH_IDX_beta[:10]-1)
# time_labels = cTROUGH_IDX_beta-1
# print(time_labels.shape)
# labels = bl1Pres[cTROUGH_IDX_beta-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# idx_run = np.where(bl1Pres[cTROUGH_IDX_beta-1]==1)
# print(idx_run[0][:10])


# # In[13]:


# # phate.plot.rotate_scatter3d(motif_beta_3D[idx_run[0],:], c=time_labels[idx_run[0],0], filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()
# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[14]:


# # # import scipy as sio
# # # fullpath = os.path.join(path,'labels_gamma_random')
# # torch.save(labels, os.path.join(path,'labels_gamma_random.pt'))
# # # sio.io.savemat('labels_gamma_random.mat', {'labels': labels})

# # torch.save(motif_beta_3D, os.path.join(path,'motif_gamma_3D_random.pt'))
# # # sio.io.savemat('motif_gamma_3D_random.mat', {'motif_gamma_3D': motif_beta_3D})


# # In[15]:




# # labels = torch.load('labels_gamma_random.pt')
# # motif_gamma_3D = torch.load('motif_gamma_3D_random.pt')
# # # motif_gamma_1D = torch.load('motif_beta_1D.pt')
# # motif_gamma = gamma_vars['cMOTIF_zscore'][0][channel]

# # from sklearn.decomposition import PCA
# # from scipy import stats
# # from sklearn import cluster
# # pca = PCA(n_components=9)


# # # print(motif_gamma_3D.shape)
# # # print(motif_gamma_1D.shape)

# # motif_gamma_run = pca.fit_transform(motif_gamma[idx_run[0],:])
# # print(motif_gamma_run.shape)

# # # phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run[0],:], c=motif_gamma_run[:,5], filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=10, alpha=0.1)
# # # plt.close()
# # # Image(filename=path + "/meso_" + 'experiment' + ".gif")

# # # print(motif_gamma_3D[idx_run[0],:].shape)
# # # kde = stats.gaussian_kde(motif_gamma_3D[idx_run[0],:])
# # # density = kde(motif_gamma_3D[idx_run[0],:])

# # # fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
# # # x, y, z = motif_gamma_3D[idx_run[0],:]
# # # ax.scatter(x, y, z, c=density)
# # # plt.show()

# # # spectral = cluster.SpectralClustering(
# # #         n_clusters=3, eigen_solver='arpack',
# # #         affinity="nearest_neighbors")


# # In[14]:


# motif_gamma_3D = torch.load(os.path.join(path,'motif_gamma_3D.pt'))


# # In[17]:


# from sklearn.decomposition import PCA
# from scipy import stats
# from sklearn import cluster


# # In[18]:


# # clustering = spectral.fit(motif_gamma_3D[idx_run[0],:])


# # In[75]:


# # clustering = cluster.SpectralClustering(n_clusters=20,
# #         eigen_solver='arpack',
# #         affinity="nearest_neighbors").fit_predict(motif_gamma_3D)


# # In[74]:



# # channel = 1 # 0-> beta; 1 -> gamma
# # del bl1Pres
# # bl1Pres = gamma_vars['bl1Run'][0]
# # cTROUGH_IDX = gamma_vars['cTROUGH_IDX'][0][channel]
# # # print(cTROUGH_IDX[:10]-1)
# # time_labels = cTROUGH_IDX-1
# # # print(time_labels.shape)
# # labels = bl1Pres[cTROUGH_IDX-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# # running=1
# # idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==running) #1-> Running, 0-> Not running
# # idx_run = idx_run[0]
# # # print(idx_run.type)
# # print(idx_run.shape)

# # print(clustering.shape)
# # # torch.save(clustering, os.path.join(path,'clustering.pt'))
# # phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run,:], c=clustering[idx_run], filename = path + "/clusters_" + 'experiment' + ".gif", figsize=(10, 8), cmap='tab20', s=10, alpha=0.1)
# # plt.close()
# # Image(filename=path + "/clusters_" + 'experiment' + ".gif")


# # In[18]:


# # Do computations per cluster
# # p = 0.0373 (the overall ratio) and standard deviation is sqrt(p*(1-p)/n) where n is the number of points in the cluster.
# channel = 1 # 0-> beta; 1 -> gamma
# del bl1Pres
# bl1Pres = gamma_vars['bl1Run'][0]
# cTROUGH_IDX = gamma_vars['cTROUGH_IDX'][0][channel]
# # print(cTROUGH_IDX[:10]-1)
# time_labels = cTROUGH_IDX-1
# # print(time_labels.shape)
# labels = bl1Pres[cTROUGH_IDX-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# running=1

# clustering = torch.load(os.path.join(path,'clustering.pt'))
# idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==1) # idx to running points
# idx_run = idx_run[0]
# sem_array = np.zeros_like(clustering, dtype=float)
# p_array = np.zeros_like(clustering, dtype=float)
# p_all = len(idx_run)/len(clustering)
# pval_array = np.zeros_like(clustering, dtype=float)
# print("p_all: {}".format(p_all))
# for i in range(0,20):
#     idx_cluster = np.where(clustering == i)
#     n = len(idx_cluster[0])
# #     print(n)
# #     print(idx_cluster[0].shape)
#     running_in_cluster = np.intersect1d(idx_cluster, idx_run)
#     p = len(running_in_cluster)/n
# #     print(p)
#     sem = np.sqrt(p_all*(1-p_all)/n)
# #     print(sem)
#     sem_array[idx_cluster[0]] = sem
#     p_array[idx_cluster[0]] = p
#     pvalue = (1 - stats.norm.cdf(np.abs(p - p_all)/sem)) * 2
# #     p_value = (1 - stats.norm.cdf(np.abs(p - p_all)/(np.sqrt(p_all*(1-p_all)/n)))) * 2 # 2 * (1 - normcdf(abs(p_cluster - p_all)/sqrt(p_all * (1 - p_all) / n_cluster) ) )
#     if pvalue==0:
#         pval_array[idx_cluster[0]] = 100
#     else:
#         pval_array[idx_cluster[0]] = -np.log(pvalue)
# #     print(-np.log(pvalue))
# #     if i > 2:
# #         break
# #     print('idx_cluster: {}'.format(idx_cluster[0][:400]))
# #     print('sem_array: {}'.format(sem_array[:400]))
# #     print('idx_run: {}'.format(idx_run[:400]))
# #     print(np.intersect1d(idx_cluster, idx_run))
#     print("cluster {}: {} points, p = {}, se = {}, pvalue= {}".format(i,n,p,sem,pvalue))
# #     break


# # In[78]:


# # # phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run,:], c=pval_array[idx_run], filename = path + "/clusters_pval_" + 'experiment' + ".gif", figsize=(10, 8), cmap='jet', vmax=100, s=10, alpha=0.1)
# # # plt.close()
# # # Image(filename=path + "/clusters_pval_" + 'experiment' + ".gif")

# # phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run,:], c=pval_array[idx_run], filename = path + "/clusters_pval_" + 'experiment' + ".gif", figsize=(10, 8), 
# #                             legend_title='-log(Pvalue)',title='Real data',cmap='jet',vmax=100, s=10, alpha=0.1)

# # # phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run,:], c=p_array[idx_run], filename = path + "/clusters_p_" + 'experiment' + ".gif", figsize=(10, 8), 
# # #                             legend_title='Ratio (p)',title='Real data',cmap='jet', s=10, alpha=0.1)

# # plt.close()
# # # Image(filename=path + "/clusters_p_" + 'experiment' + ".gif")
# # Image(filename=path + "/clusters_pval_" + 'experiment' + ".gif")


# # In[2]:


# # torch.save(motif_beta_3D, os.path.join(path,'motif_gamma_3D_random.pt'))
# path = '/home/ahf38/Documents/gamma_bouts/data'
# file = 'hilbert_vars.mat'
# fullpath = os.path.join(path,file)
# gamma_vars = sio.loadmat(fullpath)
# motif_gamma_3D = torch.load(os.path.join(path,'motif_gamma_3D.pt'))

# print('motif_gamma_3D.shape: {}'.format(motif_gamma_3D.shape))
# channel = 1 # 0-> beta; 1 -> gamma
# # del bl1Pres
# bl1Pres = gamma_vars['bl1Run'][0]
# cTROUGH_IDX = gamma_vars['cTROUGH_IDX'][0][channel]
# # print(cTROUGH_IDX[:10]-1)
# time_labels = cTROUGH_IDX-1
# # print(time_labels.shape)
# labels = bl1Pres[cTROUGH_IDX-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# running=1
# idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==running) #1-> Running, 0-> Not running
# idx_run = idx_run[0]
# # print(idx_run.type)
# print(idx_run.shape)

# print('motif_gamma_3D[idx_run,:].shape: {}'.format(motif_gamma_3D[idx_run,:].shape))
# # phate.plot.rotate_scatter3d(motif_gamma_3D_random[idx_run,:], c=labels[idx_run], filename = path + "/clusters_" + 'experiment' + ".gif", figsize=(10, 8), cmap='tab20', s=10, alpha=0.1)
# # plt.close()
# # Image(filename=path + "/clusters_" + 'experiment' + ".gif")


# # In[3]:


# # from sklearn.cluster import KMeans 
# # #We are going to use k-means for this clustering now
# # clusters_KMeans = KMeans(n_clusters = 65, max_iter = 500).fit_predict(motif_gamma_3D)


# # In[4]:


# phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run,:], c=clusters_KMeans[idx_run], filename = path + "/clusters__kmeans" + 'experiment' + ".gif", figsize=(10, 8), cmap='jet', s=10, alpha=0.1)
# plt.close()
# Image(filename=path + "/clusters__kmeans" + 'experiment' + ".gif")


# # In[6]:


# # torch.save(clusters_KMeans, os.path.join(path,'clusters_KMeans_random_data.pt'))
# clusters_KMeans = torch.load(os.path.join(path,'clusters_KMeans_data.pt'))


# # In[71]:


# # Do computations per cluster
# # p = 0.0373 (the overall ratio) and standard deviation is sqrt(p*(1-p)/n) where n is the number of points in the cluster.
# from scipy import stats
# clustering = clusters_KMeans[:]
# print(len(np.unique(clustering)))
# channel = 1 # 0-> beta; 1 -> gamma
# # del bl1Pres
# bl1Pres = gamma_vars['bl1Run'][0]
# cTROUGH_IDX = gamma_vars['cTROUGH_IDX'][0][channel]
# # print(cTROUGH_IDX[:10]-1)
# time_labels = cTROUGH_IDX-1
# # print(time_labels.shape)
# labels = bl1Pres[cTROUGH_IDX-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# running=1

# idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==1) # idx to running points
# idx_run = idx_run[0]
# sem_array = np.zeros_like(clustering, dtype=float)
# p_array = np.zeros_like(clustering, dtype=float)
# p_all = len(idx_run)/len(clustering)
# pval_array = np.zeros_like(clustering, dtype=float)
# color_cluster = np.zeros_like(clustering, dtype=float)
# centroids_and_stats = []
# list_target_clusters = [17,29,3,23]
# print("p_all: {}".format(p_all))
# for i in range(0,len(np.unique(clustering))):
#     idx_cluster = np.where(clustering == i)
#     n = len(idx_cluster[0])
# #     print(n)
# #     print(idx_cluster[0].shape)
#     running_in_cluster = np.intersect1d(idx_cluster, idx_run)
#     centroid_cluster = motif_gamma_3D[running_in_cluster,:].mean(axis=0)
#     p = len(running_in_cluster)/n
# #     print(p)
#     sem = np.sqrt(p_all*(1-p_all)/n)
# #     print(sem)
#     sem_array[idx_cluster[0]] = sem
#     p_array[idx_cluster[0]] = p
#     if i in list_target_clusters:
#         color_cluster[idx_cluster[0]]=1
#     pvalue = (1 - stats.norm.cdf(np.abs(p - p_all)/sem)) * 2
# #     p_value = (1 - stats.norm.cdf(np.abs(p - p_all)/(np.sqrt(p_all*(1-p_all)/n)))) * 2 # 2 * (1 - normcdf(abs(p_cluster - p_all)/sqrt(p_all * (1 - p_all) / n_cluster) ) )
#     if pvalue==0:
#         pval_array[idx_cluster[0]] = 100
#         print("cluster {}: {} points, p = {:.6f}, se = {:.6f}, pvalue= {:.6f}, cent = {}".format(i,n,p,sem,pvalue, centroid_cluster))
#         centroids_and_stats.append([np.array(centroid_cluster)[0], np.array(centroid_cluster)[1], np.array(centroid_cluster)[2], p, 100])
# #         pvalue=100
#     else:
#         pval_array[idx_cluster[0]] = -np.log(pvalue)
#         print("cluster {}: {} points, p = {:.6f}, se = {:.6f}, pvalue= {:.6f}, cent = {}".format(i,n,p,sem,pvalue, centroid_cluster))
#         centroids_and_stats.append([np.array(centroid_cluster)[0], np.array(centroid_cluster)[1], np.array(centroid_cluster)[2], p, -np.log(pvalue)])
#     # Append all 
# #     centroids_and_stats.append([np.array(centroid_cluster)[0], np.array(centroid_cluster)[1], np.array(centroid_cluster)[2], p, -np.log(pvalue)])
# #     print(-np.log(pvalue))
# #     if i > 2:
# #         break
# #     print('idx_cluster: {}'.format(idx_cluster[0][:400]))
# #     print('sem_array: {}'.format(sem_array[:400]))
# #     print('idx_run: {}'.format(idx_run[:400]))
# #     print(np.intersect1d(idx_cluster, idx_run))
# #     print("cluster {}: {} points, p = {:.6f}, se = {:.6f}, pvalue= {:.6f}, cent = {}".format(i,n,p,sem,-np.log(pvalue), centroid_cluster))
# #     break


# # In[63]:


# phate.plot.rotate_scatter3d(motif_gamma_3D[idx_run,:], c=color_cluster[idx_run], filename = path + "/clusters_kmeans_hierarchical_" + 'experiment' + ".gif", figsize=(10, 8), 
#                             legend_title='-log(Pvalue)',title='Real data',cmap='jet',vmax=100, s=10, alpha=0.1)

# # phate.plot.rotate_scatter3d(motif_gamma_3D_random[idx_run,:], c=p_array[idx_run], filename = path + "/clusters_kmeans_p_" + 'experiment' + ".gif", figsize=(10, 8), 
# #                             legend_title='Ratio (p)',title='Random data',cmap='jet', s=10, alpha=0.1)

# plt.close()
# # Image(filename=path + "/clusters_kmeans_p_" + 'experiment' + ".gif")
# Image(filename=path + "/clusters_kmeans_hierarchical_" + 'experiment' + ".gif")


# # In[73]:


# # Use hierarchical clusetring to know which groups should be merged together
# import scipy.cluster.hierarchy as shc
# from sklearn.preprocessing import normalize
# # centroids_and_stats = normalize(centroids_and_stats)
# # Data to be used: p (ratio), -log(P_val) and centroids of the clusters
# centroids_and_stats = np.array(centroids_and_stats)
# centroids_and_stats = centroids_and_stats / centroids_and_stats.max(axis=0)
# print(centroids_and_stats.shape)
# print(centroids_and_stats)
# plt.figure(figsize=(15, 15))  
# plt.title("Dendrograms")  
# dend = shc.dendrogram(shc.linkage(centroids_and_stats, method='ward'), leaf_font_size=15)


# # In[12]:


# #Redo embedding with UMAP
# path = '/home/ahf38/Documents/gamma_bouts/data'
# file = 'hilbert_vars.mat'
# fullpath = os.path.join(path,file)
# gamma_vars = sio.loadmat(fullpath)
# # motif_gamma = torch.load(os.path.join(path,'motif_gamma_3D.pt'))

# channel = 1 # 0-> beta; 1 -> gamma
# print('Channel: {}'.format(channel))
# motif_gamma = gamma_vars['cMOTIF_zscore'][0][channel]

# reducer = umap.UMAP(n_components=3, n_neighbors=5)
# data_umap = reducer.fit_transform(motif_gamma)

# # gamma=0
# # phate_operator = phate.PHATE(n_components=2,gamma=gamma, n_jobs=-1)
# # motif_beta_2D = phate_operator.fit_transform(motif_beta)


# # In[4]:


# bl1Pres = gamma_vars['bl1Run'][0]
# cTROUGH_IDX = gamma_vars['cTROUGH_IDX'][0][channel]
# # print(cTROUGH_IDX[:10]-1)
# time_labels = cTROUGH_IDX-1
# # print(time_labels.shape)
# labels = bl1Pres[cTROUGH_IDX-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# running=1

# idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==1) # idx to running points
# idx_run = idx_run[0]


# # In[13]:


# phate.plot.rotate_scatter3d(data_umap[:,:], c=labels, filename = path + "/embed_umap_neighbors5_" + 'experiment' + ".gif", figsize=(10, 8), 
#                             legend_title='Running',title='Real data',cmap='jet', s=10, alpha=0.1)

# plt.close()
# # Image(filename=path + "/clusters_kmeans_p_" + 'experiment' + ".gif")
# Image(filename=path + "/embed_umap_neighbors5_" + 'experiment' + ".gif")


# # In[14]:


# phate.plot.rotate_scatter3d(data_umap[idx_run,:], c=labels[idx_run,:], filename = path + "/embed_umap_neighbors5_run_" + 'experiment' + ".gif", figsize=(10, 8), 
#                             legend_title='Running',title='Real data',cmap='jet', s=10, alpha=0.1)

# plt.close()
# # Image(filename=path + "/clusters_kmeans_p_" + 'experiment' + ".gif")
# Image(filename=path + "/embed_umap_neighbors5_run_" + 'experiment' + ".gif")


# # In[25]:


# # motif_gamma_3D = torch.load(os.path.join(path,'motif_gamma_3D.pt'))
# # labels = torch.load('labels_gamma_random.pt')
# # labels = torch.load(os.path.join(path,'labels_gamma_random.pt'))
# # # motif_gamma_3D = torch.load('motif_gamma_3D_random.pt')
# # # motif_gamma_1D = torch.load('motif_beta_1D.pt')
# # motif_gamma = gamma_vars['cMOTIF_zscore'][0][channel]


# # In[79]:


# # channel = 1 # 0-> beta; 1 -> gamma
# # del bl1Pres
# # bl1Pres = gamma_vars['bl1Run'][0]
# # cTROUGH_IDX = gamma_vars['cTROUGH_IDX'][0][channel]
# # # print(cTROUGH_IDX[:10]-1)
# # time_labels = cTROUGH_IDX-1
# # # print(time_labels.shape)
# # labels = bl1Pres[cTROUGH_IDX-1]#Note: the idxs in cTROUGH_IDX are in matlab (starts in 1). So subtract 1 from all the indexes.
# # running=1
# # idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==running) #1-> Running, 0-> Not running
# # idx_run = idx_run[0]
# # # print(idx_run.type)
# # print(idx_run[:10])

# # # # For original dimension
# # # kde = stats.gaussian_kde(motif_gamma[idx_run[0],:].T)
# # # density = kde(motif_gamma[idx_run[0],:].T)

# # # # Dummy points for control
# # # num_points_in_running = idx_run.shape
# # # idx_run = np.random.choice(motif_gamma_3D.shape[0], size=num_points_in_running, replace=False)
# # # # random_rows = an_array[random_indices, :]
# # # print(num_points_in_running)
# # # print(idx_run.shape)
# # # print(idx_run[:10])

# # # # For embed 
# # band=0.001
# # # kde = stats.gaussian_kde(motif_gamma_3D[idx_run,:].T, bw_method=band)
# # # density = kde(motif_gamma_3D[idx_run,:].T)

# # # instantiate and fit the KDE model
# # kde = KernelDensity(bandwidth=band, kernel='gaussian')
# # kde.fit(motif_gamma_3D[idx_run,:])

# # # score_samples returns the log of the probability density
# # logprob = kde.score_samples(motif_gamma_3D[idx_run,:])
# # print(logprob.shape)
# # # density = np.log(density[:]*10**6)
# # # Plot scatter with mayavi
# # # figure = mlab.figure('DensityPlot')
# # # pts = mlab.points3d(x, y, z, density, scale_mode='none', scale_factor=0.07)
# # # mlab.axes()
# # # mlab.show()
# # # tmp = ((density - density.min())/(density.max()-density.min()))*2-1
# # # print('max:{} and min:{}'.format(tmp.max(),tmp.min()))
# # tmp1 = motif_gamma_3D[idx_run,:]
# # # idx = np.where(tmp>-10)
# # # tmp = np.log(tmp[idx])
# # # tmp1 = tmp1[idx]

# # phate.plot.rotate_scatter3d(tmp1, c=logprob, filename = path + '/real_data_' + 'band='+str(band) + ".gif", 
# #                             title='Running: '+str(running)+', band='+str(band), legend_title='PDF (LogProb)', 
# #                             figsize=(10, 8), cmap='jet', s=10, alpha=0.1)
# # plt.close()
# # Image(filename=path + '/real_data_' + 'band='+str(band) + ".gif")


# # In[14]:


# # import scipy as sio
# # # # fullpath = os.path.join(path,'labels_gamma_random')
# # # torch.save(density, os.path.join(path,'KDE_labels_sit_random.pt'))
# # # sio.io.savemat('labels_gamma_random.mat', {'labels': labels})

# # # torch.save(motif_beta_3D, os.path.join(path,'motif_gamma_3D_random.pt'))
# # sio.io.savemat(os.path.join(path, 'motif_gamma_3D.mat'), {'motif_gamma_3D': motif_gamma_3D})
# # sio.io.savemat(os.path.join(path, 'clusters_KMeans.mat'), {'clusters_KMeans': clusters_KMeans})
# # # len(density)
# # # print(tmp1.shape)


# # In[ ]:


# # tmp1 = motif_gamma_3D[idx_run[0],:]
# # # idx = np.where(tmp>-10)
# # # tmp = np.log(tmp[idx])
# # # tmp1 = tmp1[idx]

# # phate.plot.rotate_scatter3d(tmp1, c=density, filename = path + "/meso_" + 'experiment' + ".gif", 
# #                             title='Running='+str(running)+', band='+str(band), legend_title='KDE', 
# #                             figsize=(10, 8), cmap='jet', s=10, alpha=0.1)
# # plt.close()
# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[38]:


# # # Save the KDE for running and not running analyzed together
# torch.save(logprob, os.path.join(path,'logprob_run_band_'+str(band)+'.pt'))
# torch.save(motif_gamma_3D[idx_run,:], os.path.join(path,'motif_gamma_3D_run_band_'+str(band)+'.pt'))

# # torch.save(logprob, os.path.join(path,'logprob_run.pt'))
# # torch.save(motif_gamma_3D[idx_run,:], os.path.join(path,'motif_gamma_3D_run.pt'))


# # In[39]:


# # Calculate the ratio between the density of running vs not running by comparing the 
# # KDE of runnning (all points) to the KDE of not running (closest points to the KDE of running)

# # motif_gamma_3D_random = torch.load('motif_gamma_3D_random.pt')
# # # KDE_labels_sit = density[idx_run[0],:]
# # KDE_labels_random = torch.load('KDE_labels_random.pt')

# running=1
# idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==running)
# # idx_run = idx_run[0]
# motif_gamma_3D_run = torch.load(os.path.join(path, 'motif_gamma_3D_run_band_0.001.pt'))
# # KDE_labels_running = density[idx_run[0],:]
# # KDE_labels_running = torch.load('KDE_labels_running.pt')
# logprob_run = torch.load(os.path.join(path,'logprob_run_band_0.001.pt'))
# # KDE_labels_running = density 

# motif_gamma_3D_sit = torch.load(os.path.join(path, 'motif_gamma_3D_sit_band_0.001.pt'))
# running=0
# idx_run = np.where(bl1Pres[cTROUGH_IDX-1]==running)
# # idx_run = idx_run[0]
# # motif_gamma_sit = motif_gamma_3D[idx_run,:]

# # KDE_labels_sit = density[idx_run[0],:]
# # KDE_labels_sit = torch.load('KDE_labels_sit.pt')
# logprob_sit = torch.load(os.path.join(path,'logprob_sit_band_0.001.pt'))

# # For each point in running, find the closest point in sitting
# from scipy import spatial
# motif_gamma_sit_tree = motif_gamma_3D_sit
# tree = spatial.KDTree(motif_gamma_sit_tree)

# logprob_from_sit_to_run = np.zeros_like(logprob_run)
# for i in range(0,motif_gamma_3D_run.shape[0]):
# # for i in range(0,1):
# #     tree.query([(21,21)])
# #     print(motif_gamma_3D_random[i,:])
#     distance , idx_in_sit = tree.query(motif_gamma_3D_run[i,:])
# #     print('distance: {}'.format(distance))
# #     print('idx_in_sit: {}'.format(idx_in_sit))
# # #     print(tree.query(motif_gamma_running[0,:]))
# #     print(KDE_labels_running[i])
# #     print(KDE_labels_sit[idx_in_sit])
#     logprob_from_sit_to_run[i] = logprob_sit[idx_in_sit]

# logprob_run_sit_ratio = np.divide(logprob_run,logprob_from_sit_to_run)


# # In[80]:


# # phate.plot.rotate_scatter3d(motif_gamma_3D_run, c=logprob_run_sit_ratio, filename = path + "/meso_" + 'experiment' + ".gif", 
# #                             title='Ratio Random / Sitting, band='+str(band), legend_title='KDE ratio', 
# #                             figsize=(10, 8), cmap='jet', s=10, alpha=0.1)
# # plt.close()
# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[41]:


# print(KDE_running_sit_ratio[:100])


# # In[37]:


# np.exp(10)


# # In[28]:


# # Perform Kernel regression for running and not running
# # torch.save(density, 'KDE_labels_running.pt')
# import statsmodels.nonparametric.api as nparam

# KDE_labels_running = torch.load('KDE_labels_running.pt')
# print(KDE_labels_running.shape)
# print(motif_gamma_3D[idx_run[0],:].shape)
# print('Labels loaded')
# model = nparam.KernelReg(KDE_labels_running, motif_gamma_3D[idx_run[0],:], 
#                          var_type = ['c', 'c', 'c'], bw='cv_ls') 


# # In[ ]:


# mean, mfx = model.fit()
# ax = fig.subplots()
# ax.plot(f.x, mean, color='r', lw=2, label='est. mean')
# ax.legend(loc='upper left')
# res.append((model, mean, mfx))


# # In[22]:


# # density = kde(motif_gamma[idx_run[0],:].T)

# # labels = torch.load('labels.pt')
# # motif_gamma_3D = torch.load('motif_gamma_3D.pt')
# # motif_gamma_1D = torch.load('motif_beta_1D.pt')
# # motif_gamma = gamma_vars['cMOTIF_zscore'][0][channel]

# # from scipy import stats
# # kde = stats.gaussian_kde(motif_gamma[idx_run[0],:].T)
# # density = kde(motif_gamma[idx_run[0],:].T)
# # print(density.max())
# # density.min()
# tmp = ((density - density.min())/(density.max()-density.min()))*2-1
# # tmp = density[:]
# # tmp = np.sqrt(density)/density.std()
# # idx = np.where(tmp<0.5)
# # print(idx[0])
# # tmp = tmp[np.where(tmp<0.01)]
# fig, ax = plt.subplots()
# ax.hist(tmp, 1000)
# # ax.set_xlim([0,10])
# ax.set_ylim([0,100])


# # In[14]:


# print(density.sh)


# # In[ ]:


# channel = 1 # 0-> beta; 1 -> gamma
# print('Channel: {}'.format(channel))
# motif_beta = gamma_vars['cMOTIF_zscore'][0][channel]
# phate_operator = phate.PHATE(n_components=1,gamma=0, n_jobs=-1)
# motif_beta_1D = phate_operator.fit_transform(motif_beta)


# # In[ ]:


# torch.save(motif_beta_1D, 'motif_beta_1D.pt')
# sio.io.savemat('motif_gamma_3D.mat', {'motif_gamma_3D': motif_beta_3D})


# # In[31]:


# # reducer = umap.UMAP()
# # data_umap = reducer.fit_transform(motif_gamma)

# # fig, ax = plt.subplots(figsize=(10, 8), dpi= 80, facecolor='w', edgecolor='k')
# # try:
# #     im = ax.scatter(data_umap[:,0], data_umap[:,1], c=range(offset,n_samples), s=20)
# # except:
# #     im = ax.scatter(data_umap[:,0], data_umap[:,1], c=range(0,len(data_umap)), s=20)


# # In[16]:


# # indices of exclusive conditions
# # %matplotlib inline
# bl1RunOnly = gamma_vars['bl1RunOnly']
# bl1PresOnly = gamma_vars['bl1PresOnly']
# bl1Quiet = gamma_vars['bl1Quiet']

# fig, (ax,ax2) = plt.subplots(2,1,figsize=(10, 15), dpi= 80, facecolor='w', edgecolor='k')
# ax.plot(bl1RunOnly[0,:]) #There is running / whisking but not vis stim
# ax.plot(bl1PresOnly[0,:]) #There is vis stim, but not running 
# ax.plot(bl1Quiet[0,:]) #Neither stim or action
# ax.legend(['bl1RunOnly', 'bl1PresOnly', 'bl1Quiet'])
# ax.set_title('Running / whisking')
# ax.set_ylabel('Stimuli')
# ax.set_xlabel('Frame #')
# Fs = gamma_vars['inSampleRate'][0][0]  # 1kHz

# #Combine al inputs
# a = bl1RunOnly[0,:]*2
# b = bl1PresOnly[0,:] * 3
# c = bl1Quiet[0,:]
# print(len(bl1RunOnly[0,:]))
# time = range(0,len(bl1RunOnly[0,:]))/Fs

# all_stim = a+b+c
# ax2.plot(time, all_stim)


# # In[10]:


# # Compute fourier transform for each channel 
# NFFT = 200 #Number of points 
# noverlap = 4
# Fs = gamma_vars['inSampleRate'][0][0]  # 1kHz
# print('Fs: {}'.format(Fs))

# offset=0
# # video = video.reshape((video.shape[0],256, 256))
# # print(db2LFP[1,:].shape)
# fig, ax = plt.subplots(figsize=(10, 5))

# #To get dimensions
# Pxx, freqs, bins, im = ax.specgram(db2LFP[1,:], NFFT=NFFT, Fs=Fs, noverlap=noverlap, mode='phase') 
# plt.colorbar(im);
# # We can see gamma around 60Hz there

# P_phase = np.zeros((Pxx.shape[0],Pxx.shape[1],db2LFP.shape[0]))
# P_power = np.zeros((Pxx.shape[0],Pxx.shape[1],db2LFP.shape[0]))
# print(P_power.shape)

# for i in range(0, db2LFP.shape[0]):
#     print(i)
#     Pxx, freqs, bins, im = ax.specgram(db2LFP[i,:], NFFT=NFFT, Fs=Fs, noverlap=noverlap, mode='phase', scale_by_freq=True)
#     P_phase[:,:,i] = Pxx
#     Pxx, freqs, bins, im = ax.specgram(db2LFP[i,:], NFFT=NFFT, Fs=Fs, noverlap=noverlap, mode='magnitude', scale_by_freq=True)
#     P_power[:,:,i] = Pxx


# # In[11]:


# # fig, ax = plt.subplots(figsize=(10, 5))
# # Pxx, freqs, bins, im = ax.specgram(db2LFP[1,:], NFFT=NFFT, Fs=Fs, noverlap=noverlap, mode='phase') 
# # # plt.imshow(Pxx[:,:])
# # ax.set_aspect('auto')
# # plt.colorbar(im);
# # # plt.clim(0, 0.3);


# # In[12]:


# # print(freqs.shape)

# # path_freq_decomp = os.path.join(path_file,experiment,"output/freq_decomp")
# # np.savez(path_freq_decomp, P=P, freqs=freqs, bins=bins)


# # In[13]:


# print('bins.shape: {}'.format(bins.shape))
# print('bins: {}'.format(bins))
# print('freqs: {}'.format(freqs))
# print('freqs.shape: {}'.format(freqs.shape))
# # print('P_phase: {}'.format(P_phase))
# print('P_phase.shape: {}'.format(P_phase.shape))
# cutoff_freq = range(2,17)
# P_phase_cutoff = P_phase[cutoff_freq,:,:]
# P_power_cutoff = P_power[cutoff_freq,:,:]
# freqs_cutoff = freqs[cutoff_freq]
# print(freqs_cutoff)


# # In[14]:


# from scipy import stats
# fig, (ax1,ax2) = plt.subplots(2,1,figsize=(18, 15))
# ax1.plot(P_power_cutoff.max(axis=(0,2)))
# ax1.plot()
# ax1.set_title('Amplitude/sqrt(Hz)', fontsize=15)
# ax2.plot(P_phase_cutoff.max(axis=(0,2)))
# ax2.set_title('Phase', fontsize = 15)
# # P_power_cutoff.max(axis=(0,2)).shape
# # ax3.plot(stats.zscore(P_power_cutoff, axis=1))
# # ax3.set_title('Zscore power', fontsize = 15)


# # In[15]:


# print(P_phase_cutoff.shape)

# P_phase_subtracted = []
# for i in range(0,P_phase_cutoff.shape[1]):
# #     print(i)
#     P_phase_subtracted.append((P_phase_cutoff[:,i,:] - P_phase_cutoff[:,i,0].reshape(P_phase_cutoff.shape[0],1)).reshape(15,1,15))
# #     if i>2:
# #         break
    
# P_phase_subtracted = np.concatenate(P_phase_subtracted, axis=1)
# print(P_phase_subtracted.shape)

# fig, ax = plt.subplots(1,1,figsize=(18, 5), dpi= 80, facecolor='w', edgecolor='k')
# # count = np.argmin(np.abs(3500-bin[idx_time]))
# count = 3 #near stim
# print(bins[count])
# print(freqs_cutoff)
# P2 = P_phase_subtracted[:,count,:] # P[freq, time, channel]
# a = P2[:,:].transpose()
# # print(a)
# im = ax.imshow(a, interpolation='none', extent=[freqs_cutoff[0]-2,freqs_cutoff[-1]+2,0,30]) #extent=[8,82,0,14]
# fig.colorbar(im, ax=ax)
# ax.set_xticks(range(int(freqs_cutoff[0]),int(freqs_cutoff[-1])+2)[::5])
# # ax.set_yticks([])
# # ax.set_xticks([])
# ax.set_xlabel('Freq', fontsize=20)
# # ax.set_xticklabels(freqs_cutoff[:])
# ax.set_ylabel('Channel', fontsize=20)
# ax.set_aspect('auto')
# ax.set_title(str(bins[count]) + " sec", fontsize=20)


# # In[16]:


# # Combine it all in one variable
# print(P_phase_subtracted.shape)
# print(freqs)
# # for i in range():
# # P_power_cutoff = np.sqrt(P_power_cutoff)
# print(P_phase_subtracted.max())
# print(P_power_cutoff.max())

# P_power_cutoff2 = P_power_cutoff.copy()

# for i in range(0,P_power_cutoff.shape[1]):
#     P_power_cutoff2[:,i,:] = np.sqrt(P_power_cutoff[:,i,:]*(freqs_cutoff.transpose()))
    
# P_phase_subtracted2 = (P_phase_subtracted - P_phase_subtracted.min())/(P_phase_subtracted.max() - P_phase_subtracted.min())
# P_power_cutoff2 = (P_power_cutoff2 - P_power_cutoff2.min())/(P_power_cutoff2.max() - P_power_cutoff2.min())

# print(P_phase_subtracted2.shape)
# print(P_power_cutoff2.shape)

# P_all = np.concatenate((P_phase_subtracted2, P_power_cutoff2), axis=2)
# P_all.shape

# print(P_phase_subtracted2.min())
# print(P_power_cutoff2.min())


# # In[17]:


# #Downsample everything to the bins rate
# print('bins: {}'.format(bins))
# print('time: {}'.format(time))

# i=0
# idx_time = np.zeros_like(bins)
# for bin_val in bins:
#     idx_time[i] = np.argmin(np.abs(bin_val-time))
#     print(i)
#     i+=1

# print(idx_time)


# # In[18]:


# # print(idx_time.squeeze().shape)
# # time_down = time[int(idx_time.squeeze())]
# # stim_down = all_stim[int(idx_time.squeeze())]


# # In[19]:


# #For all the frequencies combined
# # print(P[:,:,count-1].squeeze().shape)
# # plt.imshow(P[:,:,1].squeeze(), aspect='auto')
# # print(np.min(P[:,:,3]))
# # %matplotlib inline
# fig, ax = plt.subplots(10,10,figsize=(18, 18), dpi= 80, facecolor='w', edgecolor='k')
# # fig1, ax1 = plt.subplots(figsize=(5, 5), dpi= 80, facecolor='w', edgecolor='k')
# # print(ax)



# count=-1
# for j in range(0,10):
#     for k in range(0,10):
# #         if j==0 and k==0:
# #             ax[j,k].imshow(tmp[10,:,:].transpose())
# # #             ax[j,k].imshow(a, interpolation='none')
# #             ax[j,k].set_yticks([])
# #             ax[j,k].set_xticks([])
# #             ax[j,k].set_title("Original")
# #         else:
#         count=count+1
# #         print(j,k,count)
#         P2 = P_phase_subtracted[count,:,:]
# #         print(P2.shape)
# #         P2 = P2.reshape(P2.shape[0],64,64)
#         a = P2[:,:].transpose()
#         ax[j,k].imshow(a, interpolation='none')
#         ax[j,k].set_yticks([])
#         ax[j,k].set_xticks([])
#         ax[j,k].set_aspect('auto')
#         ax[j,k].set_title(str(freqs[count]) + " Hz")

        


# # In[ ]:


# matplotlib.rc('xtick', labelsize=20) 
# matplotlib.rc('ytick', labelsize=20) 
# fig, (ax,ax2,ax3) = plt.subplots(3,1,figsize=(18, 15), dpi= 80, facecolor='w', edgecolor='k')
# # count = np.argmin(np.abs(3500-bin[idx_time]))
# count = 18050 #near stim
# print(bins[count])
# print(freqs_cutoff)
# P2 = P_all[:,count,:] # P[freq, time, channel]
# a = P2[:,:].transpose()
# # print(a)
# im = ax.imshow(a, interpolation='none', extent=[freqs_cutoff[0]-2,freqs_cutoff[-1]+2,0,30]) #extent=[8,82,0,14]
# fig.colorbar(im, ax=ax)
# ax.set_xticks(range(int(freqs_cutoff[0]),int(freqs_cutoff[-1])+2)[::5])
# # ax.set_yticks([])
# # ax.set_xticks([])
# ax.set_xlabel('Freq', fontsize=20)
# # ax.set_xticklabels(freqs_cutoff[:])
# ax.set_ylabel('Channel', fontsize=20)
# ax.set_aspect('auto')
# ax.set_title(str(bins[count]) + " sec", fontsize=20)


# im = ax2.imshow(a[0:14,:], interpolation='none', extent=[freqs_cutoff[0]-2,freqs_cutoff[-1]+2,0,14]) #extent=[8,82,0,14]
# fig.colorbar(im, ax=ax2)
# ax2.set_xticks(range(int(freqs_cutoff[0]),int(freqs_cutoff[-1])+2)[::5])
# # ax.set_yticks([])
# # ax.set_xticks([])
# ax2.set_xlabel('Freq', fontsize=20)
# # ax2.set_xticklabels(freqs_cutoff[:])
# ax2.set_ylabel('Channel', fontsize=20)
# ax2.set_aspect('auto')
# ax2.set_title("Phase, " + str(bins[count]) + " sec", fontsize=20)

# im = ax3.imshow(a[15:30], interpolation='none',extent=[freqs_cutoff[0]-2,freqs_cutoff[-1]+2,0,14]) #extent=[8,82,0,14]
# fig.colorbar(im, ax=ax3)
# ax.set_xticks(range(int(freqs_cutoff[0]),int(freqs_cutoff[-1])+2)[::5])
# # ax.set_yticks([])
# # ax.set_xticks([])
# ax3.set_xlabel('Freq', fontsize=20)
# # ax.set_xticklabels(freqs_cutoff[:])
# ax3.set_ylabel('Channel', fontsize=20)
# ax3.set_aspect('auto')
# ax3.set_title("Amplitude" + str(bins[count]) + " sec", fontsize=20)

# fig.tight_layout(pad=3.0)


# # In[ ]:


# print(P_phase_subtracted.min())
# print(P_phase_subtracted.max())
# print(P_phase_subtracted.mean())
# print(P_phase_subtracted.std())


# # In[ ]:


# # # Use this to save the video
# # experiment = 'quentin'

# # print(P_phase_subtracted.shape)
# # matplotlib.rcParams['animation.embed_limit'] = 2**128
# # fig, ax = plt.subplots(figsize=(18, 10), dpi= 80, facecolor='w', edgecolor='k')
# # P2 = P_phase_subtracted[:,0,:] # P[freq, time, channel]
# # a = P2[:,:].transpose()
# # print(a.shape)
# # im=ax.imshow(a, cmap='jet', interpolation='none')

# # def init():
# #     im.set_data(np.random.random((5,5)))
# #     return [im]

# # def animate(i):
# #     offset=17500
# #     print(i)
# #     count = i + offset
# #     matplotlib.rc('xtick', labelsize=20) 
# #     matplotlib.rc('ytick', labelsize=20) 
# # #     fig, ax = plt.subplots(figsize=(18, 10), dpi= 80, facecolor='w', edgecolor='k')
# #     # count = np.argmin(np.abs(3500-bin[idx_time]))
# # #     count = 18050
# # #     print(bins[count])
# #     P2 = P_phase_subtracted[:,count,:] # P[freq, time, channel]
# #     a = P2[:,:].transpose()
    
# #     ax.imshow(a, interpolation='none',extent=[8,82,0,14], vmin=-30, vmax = 30) #extent=[8,82,0,14]
# #     ax.set_xticks(range(10,81)[::5])
# #     # ax.set_yticks([])
# #     # ax.set_xticks([])
# #     ax.set_xlabel('Freq', fontsize=20)
# #     ax.set_xticklabels(freqs_cutoff[:])
# #     ax.set_ylabel('Channel', fontsize=20)
# #     ax.set_aspect('auto')
# #     ax.set_title(str(bins[count]) + " sec", fontsize=20)
# # #     fig.suptitle('Time {}'.format(bins[count], fontsize=16))
# #     im.set_array(a)
# #     return[im]

# # anim = animation.FuncAnimation(fig, animate, frames=100, interval=60, blit=True)
                 
# # path_saveplot = os.path.join(path)
# # path_name = path_saveplot + "/" + experiment + 'v2'+ '_NFFT' + str(NFFT) + '_' + str(Fs) + 'fps' + '.mp4'
# # print(path_name)
# # anim.save(path_name, writer='ffmpeg', fps=3)


# # In[ ]:


# signal = db2LFP.transpose()
# print(signal.shape)


# # In[ ]:


# # # print('\nEmbedding (from 512 to 2 or 3) with PHATE...')
# # gamma=0
# # phate_operator = phate.PHATE(n_components=2,gamma=gamma, n_jobs=-1)
# # offset = 2700000 #The first 100 frames are bad
# # n_samples = 2900000
# # data_phate_2D = phate_operator.fit_transform(signal[offset:n_samples])
# # # torch.save(data_phate_tgt, 'data_phate_2D.pt')


# # In[ ]:


# # matplotlib.rc('xtick', labelsize=15) 
# # matplotlib.rc('ytick', labelsize=15) 
# # fig, ax = plt.subplots(figsize=(10, 8), dpi= 80, facecolor='w', edgecolor='k')
# # try:
# # #     im = ax.scatter(data_phate_2D[:,0], data_phate_2D[:,1], c=range(offset,n_samples), s=20) bl1PresOnly
# #     im = ax.scatter(data_phate_2D[:,0], data_phate_2D[:,1], c=bl1PresOnly[0,offset:n_samples], s=20) 
# # except:
# #     im = ax.scatter(data_phate_2D[:,0], data_phate_2D[:,1], c=bl1PresOnly[0,offset:n_samples], s=20)

# # cbar = plt.colorbar(im)
# # cbar.set_label('Frame number', rotation=270, fontsize=20, labelpad=50)
# # ax.set_title("Raw video", fontsize=15)
# # # cb.set_label("Foo", labelpad=-1)


# # In[ ]:


# # gamma=0
# # phate_operator = phate.PHATE(n_components=3,gamma=gamma, n_jobs=-1)
# # offset = 2700000 #The first 100 frames are bad
# # n_samples = 2900000
# # data_phate_3D = phate_operator.fit_transform(signal[offset:n_samples])


# # In[ ]:


# # print(data_phate_3D.shape[0])
# # try:
# # #     phate.plot.rotate_scatter3d(data_phate_3D, c=range(offset,n_samples), filename = "./meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# #     phate.plot.rotate_scatter3d(data_phate_3D, c=bl1PresOnly[0,offset:n_samples], filename = "./meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # except:
# #     phate.plot.rotate_scatter3d(data_phate_3D, c=range(0,data_phate_3D.shape[0]), filename = "./meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()
# # # im = ax.scatter(data_phate_tgt[:,0], data_phate_tgt[:,1], data_phate_tgt[:,2], c=range(offset,n_samples), s=20)
# # # cbar = plt.colorbar(im)
# # # cbar.set_label('Frame number', rotation=270, fontsize=20)


# # In[ ]:


# # Image(filename="./meso_" + 'experiment' + ".gif")


# # In[21]:


# # Flatten the phase plots and run phate again
# # P2 = zeros([n_samples-offset,])
# # for i in range(offset,n_samples):
# fig, (ax,ax2) = plt.subplots(2,1,figsize=(18, 10))
# ax.imshow(P_all[:,18050,:].transpose().squeeze())
# print(P_all.shape)

# P2_all = np.zeros([P_all.shape[1],P_all.shape[0]*P_all.shape[2]])
# print(P2_all.shape)
# for i in range(0, P_all.shape[1]):
# #     print(i)
#     P2_all[i,:] = P_all[:,i,:].reshape(1,P_all.shape[0]*P_all.shape[2])

# #Test with this is preserving the original info
# P3_all = P2_all[18050,:].reshape([15,1,30])
# print(P3_all.shape)
# ax2.imshow(P3_all.transpose().squeeze())


# # In[22]:


# # print('\nEmbedding (from 512 to 2 or 3) with PHATE...')
# gamma=0
# phate_operator = phate.PHATE(n_components=2,gamma=gamma,t=5, n_jobs=-1)
# data_phate_2D = phate_operator.fit_transform(P2_all)
# # sio.savemat('np_vector.mat', {'vect':vect})


# # In[25]:


# import scipy as sio
# torch.save(data_phate_2D, 'data_phate_2D.pt')
# sio.io.savemat('data_phate_2D.mat', {'data_phate_2D': data_phate_2D})


# # In[137]:


# matplotlib.rc('xtick', labelsize=15) 
# matplotlib.rc('ytick', labelsize=15)

# print(idx_time.astype(int))
# # print(type(idx_time.astype(int)))
# print(bl1PresOnly[0,idx_time.astype(int)])
# fig, (ax,ax2) = plt.subplots(1,2,figsize=(10, 8), dpi= 80, facecolor='w', edgecolor='k')
# # try:
# #     im = ax.scatter(data_phate_2D[:,0], data_phate_2D[:,1], c=range(offset,n_samples), s=20) bl1PresOnly
# im = ax.scatter(data_phate_2D[:,0], data_phate_2D[:,1], c=all_stim[idx_time.astype(int)], s=20) 
# # except:
# im = ax2.scatter(data_phate_2D[:,0], data_phate_2D[:,1], c=range(0,len(idx_time)), s=20)

# # cbar = plt.colorbar(im)
# # cbar.set_label('Frame number', rotation=270, fontsize=20, labelpad=50)
# # ax.set_title("Raw video", fontsize=15)
# # cb.set_label("Foo", labelpad=-1)


# # In[26]:


# gamma=0
# phate_operator = phate.PHATE(n_components=3,gamma=gamma, n_jobs=-1)
# data_phate_3D = phate_operator.fit_transform(P2_all)
# print(data_phate_3D.shape[0])


# # In[34]:


# torch.save(data_phate_3D, 'data_phate_3D.pt')
# # print(type(data_phate_3D))
# tmp = all_stim[idx_time.astype(int)]
# # print(type(tmp))
# sio.io.savemat('data_phate_3D.mat', {'data_phate_3D': data_phate_3D})
# sio.io.savemat('labels.mat',{'labels': tmp})


# # In[35]:


# # try:
# # #     phate.plot.rotate_scatter3d(data_phate_3D, c=range(offset,n_samples), filename = "./meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# #     phate.plot.rotate_scatter3d(data_phate_3D, c=all_stim[idx_time.astype(int)], filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # except:
# #     phate.plot.rotate_scatter3d(data_phate_3D, c=range(0,data_phate_3D.shape[0]), filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()


# # In[36]:


# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[37]:


# # phate.plot.rotate_scatter3d(data_phate_3D, c=range(0,len(idx_time)), filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()


# # In[38]:


# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[143]:


# #Embedding of the points that are not associated to no behavior
# action = all_stim[idx_time.astype(int)]
# print(action.shape)
# idx_action = action[:]!=0
# action = action[idx_action]
# print(action.shape)
# data_action = P2_all[idx_action,:]
# print(data_action.shape)

# gamma=0
# phate_operator = phate.PHATE(n_components=3,gamma=gamma, n_jobs=-1)
# data_phate_3D = phate_operator.fit_transform(data_action)
# print(data_phate_3D.shape[0])


# # In[39]:


# # phate.plot.rotate_scatter3d(data_phate_3D, c=range(0,len(action)), filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()
# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[40]:


# # phate.plot.rotate_scatter3d(data_phate_3D, c=action[:], filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()
# # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[42]:


# # #Embedding of the points that are not associated to no behavior
# # idx_action = action[:]!=1
# # # action = action[idx_action]
# # print(action.shape)
# # # data_action = P2_all[idx_action,:]
# # # print(data_action.shape)

# # # gamma=0
# # # phate_operator = phate.PHATE(n_components=3,gamma=gamma, n_jobs=-1)
# # # data_phate_3D = phate_operator.fit_transform(data_action)
# # # print(data_phate_3D.shape[0])
# # phate.plot.rotate_scatter3d(data_phate_3D[idx_action,:], c=action[idx_action], filename = path + "/meso_" + 'experiment' + ".gif", figsize=(10, 8), cmap='viridis', s=20)
# # plt.close()
# # # Image(filename=path + "/meso_" + 'experiment' + ".gif")


# # In[ ]:





# # In[ ]:




