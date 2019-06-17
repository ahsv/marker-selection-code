# Alexander Vargo, Septermber 23 2018
# 
# Running MAST on the zeisel14 dataset.
# with 5000 variable genes retained

##### IMPORTS
import numpy as np
import pandas as pd

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

### USING THE ZEISEL DATA SET
data = np.load("/home/ahsvargo/publicData/zeisel/zeisel-proc.npz")

# for MAST: we want to normalize the CPMs
data = Rocks(data['X'], data['y'])
data.normalize(totalexpr=10**6, log=True)

ft = FoldTester(data)
ft.loadfolds("/home/ahsvargo/publicData/zeisel/zeisel14-5folds.npz")
ft.makerocks(1)

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


# input data that has already been log transformed
# normalized so that all rows have counts = 10**6
#
# We save the whole list of 3x3 matrices.  
# We return the list of genes in pvalue order, from lowest to highest
def MASTtpmDet_for_clust(Xtranspose,y, clust):
    
    clustGroup = 1.0 * (y==clust)

    cdr = robjects.r['scale']( np.mean(1*(Xtranspose > 0), axis=0) )
    sca = mast.FromMatrix( exprsArray=Xtranspose, 
            cData=robjects.r['data.frame']( 
                wellKey=np.array( range(Xtranspose.shape[1]) ), grp=clustGroup, cdr=cdr
                )
            )
    
    fmla = robjects.Formula("~cdr + grp")
    environment = fmla.environment
    environment['grp'] = clustGroup
    environment['cdr'] = cdr
    zlmdata = robjects.r['zlm.SingleCellAssay'](fmla, sca)
    
    result = mast.lrTest(zlmdata, "grp")

    np.savez("MASTDet-clust{}-results.npz".format(clust), np.array(result))
    
    return np.argsort(np.array(result)[:,2,2])


# all differentially expressed genes for all clusters
def MASTtpmDet(Xtranspose, y=None):
    
    diffExp = []
    for clust in range(y.max()+1):
        print("Working on cluster: {}".format(clust), flush=True)

        diffExp.append(MASTtpmDet_for_clust(Xtranspose,y, clust))
        
    return diffExp


# fold testing for MASTdet
foldData = MASTtpmDet( ft.rocks[foldN].X.T, ft.rocks[foldN].y.flatten() )
foldData = np.array(foldData)


np.savez("markerList-fold{}.npz".format(foldN), foldData)
