# Alexander Vargo, Septermber 23 2018
# 
# Running edgeR on the paul15 dataset.

##### IMPORTS
import numpy as np
import pandas as pd
import scanpy.api as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()


# pictured rocks
import sys
sys.path.append('/home/ahsvargo/xvalid')
from picturedrocks import Rocks
from picturedrocks.performance import FoldTester, PerformanceReport, NearestCentroidClassifier

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

### USING THE PAUL15 DATA SET
adata = sc.datasets.paul15()

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


# This is kind of a circular method.  We take the adata object and turn it into
# a rocks object.  Even though rocks are now using AnnData for most of their
# processing. 

lookup = list(adata.obs['paul15_clusters'].cat.categories)
y = np.array([lookup.index( adata.obs['paul15_clusters'][i] ) for i in range(adata.obs['paul15_clusters'].shape[0]) ])

data = Rocks(adata.X, y)


ft = FoldTester(data)
# honestly, the loadfolds method is kind of stupid
#ft.loadfolds("/home/ahsvargo/publicData/paul/paul15-scviFolds.npz")
ft.folds = np.load("/home/ahsvargo/publicData/paul/paul15-scviFolds.npz")["folds"]
if not ft.validatefolds(): print("Something went wrong with the folds that you loaded.")
ft.makerocks(1)


def edgeRQLF_for_clust(Xtranspose,y, clust):
    
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
foldData = edgeRQLF( ft.rocks[foldN].X.T, ft.rocks[foldN].y.flatten() )

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
