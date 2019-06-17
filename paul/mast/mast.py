# Alexander Vargo, Septermber 23 2018
# 
# Running MAST on the paul15 dataset.

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


### THIS SCRIPT IS FOR MAST
mast = importr("MAST")

### USING THE PAUL15 DATA SET
adata = sc.datasets.paul15()

# for MAST: we want to normalize the CPMs
sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**6)
sc.pp.log1p(adata)

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

# this is fast since the paul dataset is small.  It won't be smart to do this on the larger data sets.
ft = FoldTester(data)
# honestly, the loadfolds method is kind of stupid
#ft.loadfolds("/home/ahsvargo/publicData/paul/paul15-scviFolds.npz")
ft.folds = np.load("/home/ahsvargo/publicData/paul/paul15-scviFolds.npz")["folds"]
if not ft.validatefolds(): print("Something went wrong with the folds that you loaded.")
ft.makerocks(1)


# input data that has already been log transformed
# normalized so that all rows have counts = 10**6
#
# We save the whole list of 3x3 matrices.  
# We return the list of genes in pvalue order, from lowest to highest
def MASTtpm_for_clust(Xtranspose,y, clust):
    
    clustGroup = 1.0 * (y==clust)

    sca = mast.FromMatrix( exprsArray=Xtranspose, 
            cData=robjects.r['data.frame']( 
                wellKey=np.array( range(Xtranspose.shape[1]) ), grp=clustGroup 
                )
            )
    
    fmla = robjects.Formula("~grp")
    environment = fmla.environment
    environment['grp'] = clustGroup
    zlmdata = robjects.r['zlm.SingleCellAssay'](fmla, sca)
    
    result = mast.lrTest(zlmdata, "grp")

    np.savez("MAST-clust{}-results.npz".format(clust), np.array(result))
    
    return np.argsort(np.array(result)[:,2,2])


# all differentially expressed genes for all clusters
def MASTtpm(Xtranspose, y=None):
    
    diffExp = []
    for clust in range(y.max()+1):
        print("Working on cluster: {}".format(clust), flush=True)

        diffExp.append(MASTtpm_for_clust(Xtranspose,y, clust))
        
    return diffExp


# fold testing for edgeRQLF
foldData = MASTtpm( ft.rocks[foldN].X.T, ft.rocks[foldN].y.flatten() )
foldData = np.array(foldData)


np.savez("markerList-fold{}.npz".format(foldN), foldData)
