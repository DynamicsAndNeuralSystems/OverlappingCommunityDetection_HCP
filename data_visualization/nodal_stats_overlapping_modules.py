import numpy as np

def module_degree_zscore_multi(W, ci, flag=0):
    '''
    Computes the within-module degree z-score, treating overlapping module nodes 
    as belonging to a "meta-module" that consists of all nodes from all modules 
    they are part of.

    Parameters
    ----------
    W : NxN np.ndarray
        Binary/weighted directed/undirected connection matrix.
    ci : NxM np.ndarray
        Binary community affiliation matrix (N nodes × M modules),
        where ci[i, j] = 1 if node i belongs to module j.
    flag : int
        Graph type:
            0: undirected graph (default)
            1: directed graph in degree
            2: directed graph out degree
            3: directed graph in and out degree

    Returns
    -------
    Z : Nx1 np.ndarray
        Within-module degree Z-score, where overlapping nodes are evaluated 
        within a "meta-module" of all regions they are connected to.
    '''
    n, m = ci.shape  # N nodes, M modules

    if flag == 2:
        W = W.copy().T
    elif flag == 3:
        W = W.copy() + W.T

    Z = np.zeros(n)  # Initialize Z-score array

    for i in range(n):  # Iterate over each node
        # Find all modules this node belongs to
        assigned_modules = np.where(ci[i, :] == 1)[0]  
        
        if len(assigned_modules) == 0:
            continue  # Skip if node is not assigned to any module
        
        # Construct the meta-module: all nodes in the same modules as node i
        meta_module_nodes = np.where(np.any(ci[:, assigned_modules], axis=1))[0]  

        if len(meta_module_nodes) > 1:  # Skip if there's only one node
            # Submatrix for the meta-module
            sub_W = W[np.ix_(meta_module_nodes, meta_module_nodes)]
            Koi = np.sum(sub_W, axis=1)  # Degree within meta-module
            
            # Compute z-score for nodes in the meta-module
            mu_K = np.mean(Koi)
            sigma_K = np.std(Koi)
            if sigma_K > 0:
                Z[i] = (Koi[np.where(meta_module_nodes == i)] - mu_K) / sigma_K

    Z[np.isnan(Z)] = 0  # Handle NaNs
    return Z

import numpy as np

def participation_coef_overlap(W, ci, degree='undirected'):
    '''
    Computes the participation coefficient for nodes in a network with overlapping modules,
    ensuring that nodes assigned to multiple modules do not artificially inflate P.

    Parameters
    ----------
    W : NxN np.ndarray
        Binary/weighted directed/undirected connection matrix.
    ci : NxM np.ndarray
        Binary community affiliation matrix (N nodes × M modules),
        where ci[i, j] = 1 if node i belongs to module j.
    degree : str
        Graph type:
            'undirected' : For undirected graphs (default).
            'in' : Uses the in-degree.
            'out' : Uses the out-degree.

    Returns
    -------
    P : Nx1 np.ndarray
        Participation coefficient.
    '''
    if degree == 'in':
        W = W.T

    n, m = ci.shape  # N nodes, M modules
    Ko = np.sum(W, axis=1)  # (out) degree per node
    Kc2 = np.zeros(n)  # community-specific degree squared

    # Compute proportional degree contribution per module
    for j in range(m):  # Iterate over modules
        module_nodes = np.where(ci[:, j] == 1)[0]  # Nodes in module j
        if len(module_nodes) == 0:
            continue
        
        # Identify neighbors that are also in module j
        sub_W = W[:, module_nodes]  # Only consider connections to this module
        Kc = np.sum(sub_W, axis=1)  # Degree of node *within* module j

        # Normalize by number of modules each node belongs to
        module_counts = np.sum(ci, axis=1)  # How many modules each node is in
        normalized_Kc = Kc / module_counts  # Distribute degree proportionally

        Kc2 += np.square(normalized_Kc)  # Sum squared contributions

    P = np.ones(n) - Kc2 / np.square(Ko)  # Compute participation coefficient
    P[np.where(Ko == 0)] = 0  # Set P=0 for isolated nodes

    return P
