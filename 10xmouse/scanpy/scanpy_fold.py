# Alexander Vargo
# 16 February 2019
#
# Running the scanpy methods on the 1M mouse data.
#
# Just changed so that we can import - no idea if this will work now

import numpy as np
import pandas as pd
import scanpy.api as sc

from scipy.stats import rankdata
import scipy.sparse as spsp

## define default fold and method method
foldN = 0
methods = ['wilcoxon', 't-test_overestim_var', 'logreg']
method = methods[1]

### READ PARAMETERS FROM FILES
ARRAY_JOB = True

if ARRAY_JOB:
    # Get the array ID and switch into the correct directory for an array job
    # TODO - probably eventually use a context manager.  
    import os
    arrID = os.getenv('PBS_ARRAYID', 1)

    path = os.getcwd()
    path = path + "/fold" + str(arrID)
    os.chdir(path)
    print("Current working directory: {}".format(os.getcwd()))

    # read in the max sparsity parameter from the working directory
    pFilename = "foldNum.dat"
    pFile = open(pFilename, 'r')
    foldN = pFile.readline().strip()
    foldN = int(foldN)

    methodN = pFile.readline().strip()
    method = methods[int(methodN)]

## needed helper function
def geneName2index(adata, name):
    geneNames = np.array(adata.var.index)
    inds = np.where(geneNames == name)[0]
    if inds.size == 0:
        print("Waring: gene name not found in adata variable index: returning 0", flush=True)
        return 0
    else:
        return inds[0]


# mark this if you are comparing only two clusters.
twoGroups = False
def scanpyMarkers(foldN, method, path=None, twoGroups=False):

    print("Finding markers for fold" + str(foldN), flush=True)
    print("Method is " + method, flush=True)

    if path == None:
        path = "/home/ahsvargo/turbo/scData/10x/scanpy/fold{}".format(foldN)

    ### LOAD DATA
    #folds = np.load("10x-5folds.npz")
    #folds = [folds["fold{}".format(i)] for i in range(5)]

    dName = path + "/1M-fold{}.h5ad".format(foldN,foldN)
    #dName = "/home/ahsvargo/turbo/scData/10x/1M-nzGenes-clusts.h5ad"
    print("Loading data from file " + dName + "...", flush=True)
    fold_adata = sc.read_h5ad(dName)
    print("done", flush=True)
    sc.logging.print_memory_usage()
    print("Number of cells in fold {}: {}".format(foldN, fold_adata.X.shape[0]), flush=True)

    ### PARAMETERS FOR rank_genes_groups
    nMarkers = 24015
    groupby = "louvain"             # louvain or graphclust

    ### RUN THE METHOD
    pvals = []
    marks = []
    clustOrder = []

    # run the method on the fold
    if twoGroups: groups = [0,1]
    else: groups = 'all' 

    print("Running method for fold {}...".format(foldN))
    sc.tl.rank_genes_groups(
        fold_adata,
        groupby=groupby,
        n_genes=nMarkers,
        method=method,
        rankby_abs=True,
        groups=groups,
        solver='saga',
        n_jobs=-1
    )
    print("done fold {}".format(foldN), flush=True)
    sc.logging.print_memory_usage()

    # initialize the saved data so that I can refer to it by index later.
    for index, name in enumerate(fold_adata.uns['rank_genes_groups']['names'].dtype.names):
        pvals.append([])
        marks.append([])
        clustOrder.append(name)
        
        pvals[index] = fold_adata.uns['rank_genes_groups']['scores'][name]
        currMarks = fold_adata.uns['rank_genes_groups']['names'][name]
        marks[index] = np.array([ geneName2index(fold_adata, gene) for gene in currMarks ])
        if (method == 'logreg' and twoGroups):
            pvals.append(pvals[0])
            marks.append(marks[0])

    pvals = np.array(pvals)
    marks = np.array(marks)

    return pvals, marks, clustOrder

def main(foldN, method, path=None, twoGroups=False):

    pvals, marks, clustOrder = scanpyMarkers(foldN, method, path, twoGroups)

    ### Save Data
    oName = "1M-fold{}-{}.npz".format(foldN, method)
    np.savez(oName, marks=marks, pvals=pvals, names=clustOrder)

if __name__ == "__main__":
    if ARRAY_JOB:
        main(foldN, method, path=path, twoGroups=twoGroups)
    else:
        main(foldN, method, path=None, twoGroups=twoGroups)
