# CBASS_L2_Embed.py;
'''
This function performs the embedding of the data specified by variable 'experiment' in the temporary directory.
The current methods available for embedding are UMAP, PHATE and PCA. Check description of the Matlab Wrapper function
for further details on the inputs: https://github.com/ahof1704/gamma_bouts/blob/master/Pipeline/CBASS_L2_Embed.m
'''
import os
import numpy as np
import umap
from sklearn.decomposition import PCA
import scipy.io as sio
from scipy import stats

#For plotting the video
import matplotlib
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--tmpPath', type=str)
parser.add_argument('--outPath', type=str)
parser.add_argument('--n_components', type=int)
parser.add_argument('--label_name', type=str)
parser.add_argument('--method', type=str)
parser.add_argument('--experiment', type=str)
args = parser.parse_args()

print('loading {}'.format(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '.mat')))
mat_vars = sio.loadmat(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '.mat'))

#Delete the temp file
# os.remove(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '.mat'))

#Perform embedding with chosen method
if args.method=='umap':
    print('Running UMAP ({}D)'.format(args.n_components))
    reducer = umap.UMAP(n_components=args.n_components, n_neighbors=5)
    emb_data = reducer.fit_transform(mat_vars['data'])
    sio.savemat(os.path.join(args.outPath,'emb_data_umap_' + args.experiment + '.mat'), {'emb_data': emb_data})

elif args.method=='pca':
    print('Running PCA ({}D)'.format(args.n_components))
    pca = PCA(n_components=args.n_components)
    emb_data = pca.fit_transform(mat_vars['data'])
    sio.savemat(os.path.join(args.outPath,'emb_data_pca_' + args.experiment + '.mat'), {'emb_data': emb_data})
    
else:
    import phate
    print('Running PHATE ({}D)'.format(args.n_components))
    phate_operator = phate.PHATE(n_components=args.n_components,gamma=0, n_jobs=-1)
    emb_data = phate_operator.fit_transform(mat_vars['data'])

    sio.savemat(os.path.join(args.outPath,'emb_data_phate_' + args.experiment + '.mat'), {'emb_data': emb_data})