import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def L(A, normalized=True, verbose=False):
    """L: compute a graph laplacian

    Args:
        A (N x N np.ndarray): Adjacency matrix of graph
        normalized (bool, optional): Normalized or combinatorial Laplacian

    Returns:
        L (N x N np.ndarray): graph Laplacian
    """
    # Compute the degree as D = diag(d(i)), d(i)=sum_{j=1}^N W_ij
    if verbose: print('A: {}'.format(A))
    D = np.sum(A, axis=0) # row sum
    if verbose: print('D: {}'.format(D))
    L = np.diag((D))-A
    if normalized:
        D_right = np.diag((D)**-0.5)
        D_left = np.diag((D)**-0.5)
        L = np.matmul(D_right, np.matmul(L,D_left)) # Normalized Markov matrix
    if verbose: print('L: {}'.format(L))
    return L
    
def SC(L, k, psi=None, nrep=5, itermax=300, show_plot=True, topK = 5, verbose=False):
    """SC: Perform spectral clustering via the Ng method
    Args:
        L (np.ndarray): Normalized graph Laplacian
        k (integer): number of clusters to compute. If k=None, estimate k based on eigengap
        nrep (int): Number of repetitions to average for final clustering
        itermax (int): Number of iterations to perform before terminating
        show_plot: Plot eigenvalues if True
        topK: shows up to K optimal number of clusters
    Returns:
        labels (N x 1 np.array): Learned cluster labels
    """
    if psi is None:
        # compute the first k elements of the Fourier basis
        e, psi = np.linalg.eigh(L) # Compute eigendecomposition
        psi_k = psi[:, :k]

    else:  # just grab the first k eigenvectors
        psi_k = psi[:, :k]

    # normalize your eigenvector rows
    psi_norm = psi_k / np.linalg.norm(psi_k, axis=1, keepdims=True)
    if verbose:
        print('psi_k: {}'.format(psi_k))
        print('psi_norm: {}'.format(psi_norm))
        print('np.sum(psi_norm[0,:]**2): {}'.format(np.sum(psi_norm[0,:]**2)))
    
    if show_plot:
        fig,ax=plt.subplots(figsize=(10,10))
        plt.title('Largest eigenvalues of input matrix')
        plt.scatter(np.arange(1, 1+len(e)), e)
        plt.xlabel("Eigenval")
        plt.ylabel("Value")
        plt.grid()
        
    # Identify the optimal number of clusters as the index corresponding
    # to the larger gap between eigen values
    inLargestGap = np.argsort(np.diff(e))[::-1][:topK]
    inNclusters = inLargestGap + 1
    if verbose: print('Optimal number of clusters: ', inNclusters)
        
    if k is None:
        k = inNclusters[0] #Get the index for the largest gap

    labels = KMeans(n_clusters=k, n_init=nrep,  max_iter=itermax).fit_predict(psi_norm)
    
    return labels, inNclusters
