# Alexander Vargo, Septermber 23 2018
# 
# Running edgeR on the paul15 dataset.

##### IMPORTS
import numpy as np
import pandas as pd
import scanpy.api as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

# rpy2
import rpy2
print(rpy2.__version__)

import rpy2.robjects as robjects

from rpy2.robjects.packages import importr

# import R's "base" package
base = importr('base')

# import R's "utils" package
utils = importr('utils')

# import rpy2's package module
import rpy2.robjects.packages as rpackages

# automatically convert numpy arrays into R data types
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects import pandas2ri
pandas2ri.activate()


### THIS SCRIPT IS FOR EDGER
edgeR = importr("edgeR")

### USING THE ZHENGFILT DATA SET
labels = "bulk"

path= "/home/ahsvargo/turbo/scData/zheng17/filtered_matrices_mex/hg19/"
foldPath = path + "zheng17-5folds.npz"

adata = sc.read_h5ad(path + "write/raw_data_5000_genes_with_clusters.h5ad")
foldFile = np.load(foldPath)
folds = [foldFile['fold' + str(i)] for i in range(5)]

# create integer cluster labels
bulk_lookup = list(np.unique(adata.obs['bulk_labels'].values))
louv_lookup = [str(a) for a in np.unique(adata.obs['louv_labels'].values)]

bulky = np.array(
        [ bulk_lookup.index( adata.obs['bulk_labels'][i]  ) 
            for i in range(adata.obs['bulk_labels'].shape[0]) ]
        ) 
louvy = np.load(path + 'zheng17_yVec_lvals.npz')['y']

if labels is 'bulk':
    yVec = bulky
elif labels is 'louvain':
    yVec = louvy


### READ PARAMETERS FROM FILES
ARRAY_JOB = True

if ARRAY_JOB:
    # Get the array ID and switch into the correct directory for an array job
    # TODO - probably eventually use a context manager.  
    import os
    arrID = os.getenv('PBS_ARRAYID', 1)

    path = os.getcwd()
    path = path + "/fold" + arrID
    os.chdir(path)
    print("Current working directory: {}".format(os.getcwd()))

    # read in the max sparsity parameter from the working directory
    pFilename = "foldNum.dat"
    pFile = open(pFilename, 'r')
    foldN = pFile.readline().strip()
    foldN = int(foldN)
else:
    foldN = 0


# get the data for the fold
N = adata.X.shape[0]
mask = np.zeros(N, dtype=bool)
mask[folds[foldN]] = True

X = np.asarray(adata.X[~mask].todense())
y = yVec[~mask]


def edgeRQLF_for_clust(Xtranspose,y, clust):
    
    print("Running for cluster {}".format(clust), flush=True)
    clustGroup = 1.0 * (y==clust)
    dge = edgeR.DGEList(counts=Xtranspose, group=clustGroup)
    dge = edgeR.calcNormFactors(dge)
    
    fmla = robjects.Formula("~grps")
    environment = fmla.environment
    environment['grps'] = clustGroup
    design = robjects.r['model.matrix'](fmla)
    
    dge = edgeR.estimateDisp(dge, design=design)
    fit = edgeR.glmQLFit(dge, design=design)
    qlf = edgeR.glmQLFTest(fit)
    
    tt = edgeR.topTags(qlf, n=Xtranspose.shape[0])
    return tt


def edgeRQLF(Xtranspose, y=None):
    
    # all differentially expressed genes for all clusters
    diffExp = []
    for clust in range(y.max()+1):
        print(clust)
        
        diffExp.append(edgeRQLF_for_clust(Xtranspose,y, clust))
        
    return diffExp


# fold testing for edgeRQLF
foldData = edgeRQLF( X.T, y.flatten() )

# write the full data out to files
for i in range(y.max()+1):
    test=pandas2ri.ri2py_dataframe(foldData[i][0])
    test.index = np.array([int(gene)-1 for gene in list(robjects.r['row.names'](foldData[i][0])) ])
    test.to_csv("fold{}-clust{}.csv".format(foldN,i))

markerList = []
for i in range(y.max()+1):
    markerList.append(
            np.array([int(gene)-1 for gene in list(robjects.r['row.names'](foldData[i][0])) ])
        )

markerList = np.array(markerList)

np.savez("markerList-fold{}.npz".format(foldN), markerList)
