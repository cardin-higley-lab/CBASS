# CBASS_L2_embedding.py;


import os
import numpy as np
import phate
#import torch
#import torch.nn as nn
import umap
from sklearn.decomposition import PCA
import scipy.io as sio
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
args = parser.parse_args()

mat_vars = sio.loadmat(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '.mat'))

#Delete the temp file
os.remove(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '.mat'))

#Perform embedding with chosen method
if args.method=='umap':
    print('Running UMAP ({}D)'.format(args.n_components))
    reducer = umap.UMAP(n_components=args.n_components, n_neighbors=5)
    emb_data = reducer.fit_transform(mat_vars['data'])
    sio.savemat(os.path.join(args.outPath,'emb_data_umap_' + args.experiment + '.mat'), {'emb_data': emb_data})

if args.method=='pca':
    print('Running PCA ({}D)'.format(args.n_components))
    pca = PCA(n_components=args.n_components)
    emb_data = pca.fit_transform(mat_vars['data'])
    sio.savemat(os.path.join(args.outPath,'emb_data_pca_' + args.experiment + '.mat'), {'emb_data': emb_data})
    
else:
    print('Running PHATE ({}D)'.format(args.n_components))
    phate_operator = phate.PHATE(n_components=args.n_components,gamma=0, n_jobs=-1)
    emb_data = phate_operator.fit_transform(mat_vars['data'])

    sio.savemat(os.path.join(args.outPath,'emb_data_phate_' + args.experiment + '.mat'), {'emb_data': emb_data})