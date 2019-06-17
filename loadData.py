# April 27 2019 

# Some functions to help with loading the various data sets that we are looking
# at.
#
# Question: will I want to load an adata object?  Rocks object?  Just the data?
#
# Probably just the data, then I can decide what I want to do with it.  Some
# data sets have multiple clusterings, however...
#
# This currently contains functions for Paul, Ziesel, Zheng, and Zhengsim
# but the zhengsim functionality hasn't been tested.
#
# Need to change the path to the data and the folds if you don't have it at 
# the same place as me (you don't)

import numpy as np
import pandas as pd
import scipy.sparse as spsp


# Load PicturedRocks for returning Rocks
import sys

# path to modified pictured rocks
sys.path.append('/home/ahsvargo/xvalid')
from picturedrocks import Rocks
from picturedrocks.performance import FoldTester, PerformanceReport, NearestCentroidClassifier

# Load scanpy for returning adata objects...
import scanpy.api as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
import anndata


## Import Data

### For Paul:

def load_paul(returnT=None):

    dataPath = "/home/ahsvargo/publicData/data/paul15/paul15.h5"
    foldPath = "/home/ahsvargo/publicData/paul15-scviFolds.npz"
    folds = np.load(foldPath)["folds"]

    adata = sc.datasets.paul15()
    sc.logging.print_memory_usage()

    lookup = list(adata.obs['paul15_clusters'].cat.categories)
    yVec = np.array([
        lookup.index( adata.obs['paul15_clusters'][i] ) 
        for i in range(adata.obs['paul15_clusters'].shape[0]) 
        ])

    if returnT is 'raw':
        return adata.X, yVec, folds

    elif returnT is 'Rocks':
        data = Rocks(adata.X, yVec)
        ft = FoldTester(data)
        ft.folds = folds
        ft.validatefolds()
        ft.makerocks(0)

        return data, ft

    elif returnT is 'adata':
        return adata, folds

    else:
        print("Return type needs to be 'raw', 'Rocks', or 'adata'")



### For zeisel:
def load_zeisel(returnT=None):

    dataPath = "/home/ahsvargo/publicData/zeisel/zeisel-proc.npz"
    foldPath = "/home/ahsvargo/publicData/zeisel/zeisel14-5folds.npz"
    data = np.load(dataPath)
    folds = np.load(foldPath)

    if returnT is 'raw':
        foldReturn = [folds['fold' + str(i)] for i in range(5)]
        return data['X'], data['y'], foldReturn

    elif returnT is 'Rocks':
        data = Rocks(data['X'], data['y'])
        ft = FoldTester(data)
        ft.loadfolds(foldPath)
        ft.makerocks(verbose=0)
        
        return data, ft

    elif returnT is 'adata':

        adata = anndata.AnnData(X = data['X'])
        adata.obs['zeisel_clusters'] = pd.Series(data['y'], dtype="category", index=adata.obs.index)
        foldReturn = [folds['fold' + str(i)] for i in range(5)]

        return adata, foldReturn

    else:
        print("Return type needs to be 'raw', 'Rocks', or 'adata'")
        


### For Zheng:
def load_zheng(returnT='None', allGenes=True, labels='bulk'):

    path= "/home/ahsvargo/turbo/scData/zheng17/filtered_matrices_mex/hg19/"
    foldPath = path + "zheng17-5folds.npz"
    folds = np.load(foldPath)

    if allGenes:
        ## Load all genes
        adata = sc.read_h5ad(path + "write/raw_data_all_nz_genes_with_clusters.h5ad")

    else:   
        ## load top 5000 genes
        adata = sc.read_h5ad(path + "write/raw_data_5000_genes_with_clusters.h5ad")

    # create integer cluster labels
    bulk_lookup = list(np.unique(adata.obs['bulk_labels'].values))
    louv_lookup = [str(a) for a in np.unique(adata.obs['louv_labels'].values)]

    bulky = np.array(
            [ bulk_lookup.index( adata.obs['bulk_labels'][i]  ) 
                for i in range(adata.obs['bulk_labels'].shape[0]) ]
            ) 
    louvy = np.load(path + 'zheng17_yVec_lvals.npz')['y']

    if labels is 'bulk':
        y = bulky
    elif labels is 'louvain':
        y = louvy
    else:
        print("Cluster labels need to be either 'bulk' or 'louvain'; returning nothing")
        return

    if returnT is 'raw':
        foldReturn = [folds['fold' + str(i)] for i in range(5)]
        return adata.X, y, foldReturn

    elif returnT is 'Rocks':
        data = Rocks(adata.X, y)
        ft = FoldTester(data)
        ft.folds = [folds['fold' + str(i)] for i in range(5)]
        ft.validatefolds()
        ft.makerocks(0)

        return data, ft

    elif returnT is 'adata':
        print("Return object contains bulk labels under obs.('bulk_labels')")
        print("and louvain labels under obs.('louv_labels')")

        foldReturn = [folds['fold' + str(i)] for i in range(5)]
        return adata, foldReturn

    else:
        print("Return type needs to be 'raw', 'Rocks', or 'adata'")



### Zheng Simulated
def makePath(test=0, allQ=True, filterQ=False):

    if allQ:
        pathPath = "allGenes"
    if filterQ and not allQ:
        partPath = "filtered"

    fname = "sparseCounts.npz"
    if allQ and filterQ:
        fname = "filtered-splat.h5ad"

    path = "/home/ahsvargo/turbo/scData/zheng17/splatter-bCells/" + partPath + "/" + str(test) + "/"
    return path, fname


def load_zhengsim(returnT='None', test=0, allQ=True, filterQ=False):

    path, fname = makePath(test=test, allQ=allQ, filterQ=filterQ)

    y = np.loadtxt(path + "y.dat", dtype='int')
    y = y - 1

    if "npz" in fname:
        X = spsp.load_npz(path + fname)
    elif "h5ad" in fname:
        adata = sc.read_h5ad(path + fname)
        X = adata.X

    foldPath = path + "bCells-5folds.npz"
    folds = np.load(foldPath)

    if returnT is 'raw':
        foldReturn = [folds['fold' + str(i)] for i in range(5)]
        return X, y, foldReturn

    elif returnT is 'Rocks':
        data = Rocks(X,y)
        ft = FoldTester(data)
        ft.loadfolds(foldPath)
        ft.makerocks(verbose=0)

        return data, ft
    
    elif returnT is 'adata':
        adata = anndata.AnnData(X = X)
        adata.obs['zhengsim_clusters'] = pd.Series(y, dtype="category", index=adata.obs.index)
        foldReturn = [folds['fold' + str(i)] for i in range(5)]

        return adata, foldReturn

    else:
        print("Return type needs to be 'raw', 'Rocks', or 'adata'")



