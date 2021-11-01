# CBASS_L2_PLotEmbed.py;
'''
This function formats and plots the embedding of the data specified by variable 'experiment' in the temporary directory.
The current method makes use of the scprep (https://scprep.readthedocs.io/en/stable/index.html), which should be automatically installed 
when replicating our conda environment (https://github.com/ahof1704/gamma_bouts/wiki/Tutorial#3-python-integration-for-embedding-and-plots).
Refer to wrapper function for further details on the inputs: https://github.com/ahof1704/gamma_bouts/blob/master/Pipeline/CBASS_L2_PlotEmbed.m
'''

import os
import numpy as np
import scprep
import matplotlib.pyplot as plt
import umap
from sklearn.decomposition import PCA
import scipy.io as sio
from scipy import stats
import matplotlib
from matplotlib import animation
import argparse
from IPython.display import HTML

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument('--tmpPath', type=str)
parser.add_argument('--outPath', type=str)
parser.add_argument('--n_components', type=int)
parser.add_argument('--method', type=str)
parser.add_argument('--experiment', type=str)
parser.add_argument('--label_name', type=str)
parser.add_argument('--add_legend', type=str2bool) 
parser.add_argument('--fontsize', type=int) 
parser.add_argument('--discrete', type=str2bool) 
parser.add_argument('--category', type=int)
parser.add_argument('--format_img', type=str)
parser.add_argument('--rotate3D', type=str2bool)
args = parser.parse_args()

matplotlib.use('Agg')

mat_vars = sio.loadmat(os.path.join(args.tmpPath,'tmp_vars_' + args.experiment + '_' + args.label_name + '.mat'))
emb_data = sio.loadmat(os.path.join(args.outPath,'emb_data_' + args.method + '_' + args.experiment + '.mat'))['emb_data']

if '20Clusters' in args.label_name:
    cmap = 'tab20'
else:
    cmap = 'viridis'

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
    
if args.n_components == 3:    
	if args.rotate3D == True:
		scprep.plot.rotate_scatter3d(emb_data, c=labels, filename = args.outPath + '/plot_embedding_' + str(args.n_components) + 'D_' + args.method + '_' + args.experiment + '_' + args.label_name + '.gif', 
		                             figsize=(10, 8), discrete = args.discrete, fontsize = args.fontsize, legend = args.add_legend, colorbar = args.add_legend, cmap=cmap, s=10, alpha=0.1)
	else:
		if args.format_img == 'eps':
			scprep.plot.scatter3d( emb_data, c=labels, filename = args.outPath + '/plot_embedding_' + str(args.n_components) + 'D_' + args.method + '_' + args.experiment + '_' + args.label_name  + '.' + args.format_img, 
				figsize=(10, 8), discrete = args.discrete, fontsize = args.fontsize, legend = args.add_legend, colorbar = args.add_legend, cmap=cmap, s=10)
		else:
			scprep.plot.scatter3d( emb_data, c=labels, filename = args.outPath + '/plot_embedding_' + str(args.n_components) + 'D_' + args.method + '_' + args.experiment + '_' + args.label_name  + '.' + args.format_img,
				dpi = 1000, discrete = args.discrete, colorbar = args.add_legend, fontsize = args.fontsize, figsize=(10, 8), cmap=cmap, s=10, alpha=0.1)
else:
	if args.format_img == 'eps':
		scprep.plot.scatter2d(emb_data, c=labels, filename = args.outPath + '/plot_embedding_' + str(args.n_components) + 'D_' + args.method + '_' + args.experiment + '_' + args.label_name + '.' + args.format_img, 
		                             figsize=(10, 8), discrete = args.discrete, fontsize = args.fontsize, legend = args.add_legend, colorbar = args.add_legend, cmap=cmap, s=10)
	else:
		scprep.plot.scatter2d(emb_data, c=labels, filename = args.outPath + '/plot_embedding_' + str(args.n_components) + 'D_' + args.method + '_' + args.experiment + '_' + args.label_name + '.' + args.format_img, 
		                             dpi = 1000, discrete = args.discrete, fontsize = args.fontsize, legend = args.add_legend, colorbar = args.add_legend, figsize=(10, 8), cmap=cmap, s=10, alpha=0.1)

plt.close('all')