# CBASS_L2_embedding.py;


import os
import logging
import time
import numpy as np
import pandas as pd
import phate
import scprep
#import torch
#import torch.nn as nn
import matplotlib.pyplot as plt
import h5py
import umap
from sklearn.decomposition import PCA
import scipy.io as sio
from sklearn.neighbors import KernelDensity
from scipy import stats

#For plotting the video
import matplotlib
from matplotlib import animation
from IPython.display import HTML, Image

#from kymatio import Scattering2D
import torch.optim
import math
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--tmpPath', type=str)
parser.add_argument('--outPath', type=str)
parser.add_argument('--n_components', type=int)
parser.add_argument('--method', type=str)
parser.add_argument('--experiment', type=str)
parser.add_argument('--category', type=int)
args = parser.parse_args()

# print('Done with imports')

# print('args: {}'.format(args))

mat_vars = sio.loadmat(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '.mat'))
emb_data = sio.loadmat(os.path.join(args.outPath,'emb_data_' + args.method + '_' + args.experiment + '.mat'))['emb_data']

labels = mat_vars['labels'][0]
# time_stamps = mat_vars['time_stamps']

# labels = bl1Run[time_stamps]
# Not sure if we will want to keep it or not
# if args.category >= 0:
#     category=args.category #1-> Running, 0-> Not running
#     print('category: {}'.format(category))
#     idx_run = np.where(labels[time_stamps]==category)
#     idx_run = idx_run[0]
# else:
#     idx_run = range(0,emb_data.shape[0]) #Just frame number
    
    
scprep.plot.rotate_scatter3d(emb_data, c=labels, filename = args.outPath + '/plot_embedding_' + args.method + '_' + args.experiment + '.gif', 
                             figsize=(10, 8), cmap='viridis', s=10, alpha=0.1)
plt.close()