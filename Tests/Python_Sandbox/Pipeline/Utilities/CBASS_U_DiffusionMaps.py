def diff_map_info(W, verbose=False):
    '''
    Construct the information necessary to easily construct diffusion map for any t

    Inputs:
        W           a numpy array of size n x n containing the affinities between points

    Outputs:

        diff_vec    a numpy array of size n x n-1 containing the n-1 nontrivial eigenvectors of Markov matrix as columns
        diff_eig    a numpy array of size n-1 containing the n-1 nontrivial eigenvalues of Markov matrix

        We assume the convention that the coordinates in the diffusion vectors are in descending order
        according to eigenvalues.
    '''
    r = np.sum(W, axis=0) # row sum
    
    D_right = np.diag((r)**-0.5)
    D_left = np.diag((r)**-0.5)
    M_s = np.matmul(D_right, np.matmul(W,D_left)) # Normalized Markov matrix
    
    eigenValues, eigenVectors = np.linalg.eigh(M_s) # Compute eigendecomposition
    if verbose: 
        print('eigenValues: {}'.format(eigenValues))
    idx = eigenValues.argsort()[::-1] # Sort and invert order of eigenvalues
    idx = idx[1:]
    if verbose:
        print('idx: {}'.format(idx))
    diff_eig = eigenValues[idx]
    if verbose:
        print('diff_eig: {}'.format(diff_eig))
    diff_vec = eigenVectors[:,idx]
    if verbose:
        print('diff_vec: {}'.format(diff_vec))

    # Compute normalization
    psi = np.matmul(D_left,diff_vec)
    diff_vec = psi/np.linalg.norm(psi)
    # return the info for diffusion maps
    return diff_vec, diff_eig

def get_diff_map(diff_vec, diff_eig, t):
    '''
    Construct a diffusion map at t from eigenvalues and eigenvectors of Markov matrix

    Inputs:
        diff_vec    a numpy array of size n x n-1 containing the n-1 nontrivial eigenvectors of Markov matrix as columns
        diff_eig    a numpy array of size n-1 containing the n-1 nontrivial eigenvalues of Markov matrix
        t           diffusion time parameter t

    Outputs:
        diff_map    a numpy array of size n x n-1, the diffusion map defined for t
    '''
    diff_map = np.matmul(diff_vec,np.diag(diff_eig)**t)
    return diff_map