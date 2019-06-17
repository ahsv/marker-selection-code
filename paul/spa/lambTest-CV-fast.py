# Alexander Vargo, 5 November 2017
#
# Edited 15 March 2018 for use with Justin's NM data and a new optimization
# function that doesn't use cvx AND the picturedrocks data structures also made
# this into an array job
#
# Edited 11 June 2018 for use with picturedRocks "updates"
# NearestCentroidClassifer doesn't normalize correctly, but we compute and
# store coeffs ahead of time.  So log normalize the rocks object before
# classifying BUT ONLY ONCE
#
# Edited 29 July 2018 for use with the fast optimization methods that save a
# lot of data and then re-use those data as much as possible.  In order to
# re-use the data, we cannot change lambFac or alpha.  So we now load in (a
# group of) alphas and lambfacs from an external file.
#
# The file with parameters is "currParams.py".  It should define the following
# parameters in terms of python variables:
# lambMesh, alphaMesh, lambs
#
# This file is not for using rank-based methods at this time.
#
# Edited 16 April 2019
# For use with scanpy datasets (paul, zeisel)
# and adding RandomForest classifier

import sys
sys.path.append('/home/ahsvargo/xvalid')

import math

import pandas as pd
import numpy as np

from picturedrocks import Rocks
from picturedrocks.performance import FoldTester, NearestCentroidClassifier

import scanpy.api as sc
sc.settings.verbosity = 3
sc.logging.print_versions()

from sklearn.ensemble import RandomForestClassifier

class RandomForest:
    def __init__(self):
        self.traindata = None
        self.RFC = RandomForestClassifier(n_estimators=100, n_jobs=-1)
        
    def train(self, data):
        self.traindata = data
        self.RFC.fit(data.X, data.y[:,0])
        
    def test(self, Xtest, sparse):
        return self.RFC.predict(Xtest)
 

### READ PARAMETERS FROM FILES
ARRAY_JOB = True

if ARRAY_JOB:
    # Get the array ID and switch into the correct directory for an array job
    # TODO - probably eventually use a context manager.  
    import os
    arrID = os.getenv('PBS_ARRAYID', 1)

    path = os.getcwd()
    path = path + "/" + arrID
    os.chdir(path)
    print("Current working directory: {}".format(os.getcwd()))

    sys.path.append(path)
else:
    arrID = "noArray"

########### Read in parameters  ####################

# Update July 29 2018:
# we load in alpha and lambfac, assuming that the file is a python file
# I guess this is kind of a security risk, but... not really
from currParams import alphaMesh, lambMesh, lambs

########### Import the data  ####################

adata = sc.datasets.paul15()

lookup = list(adata.obs['paul15_clusters'].cat.categories)
yVec = np.array(
        [lookup.index( adata.obs['paul15_clusters'][i] ) 
            for i in range(adata.obs['paul15_clusters'].shape[0]) ]
        )

# We use the scviFolds in all of the other methods.
folds = np.load("paul15-scviFolds.npz")["folds"]

y = yVec
X = adata.X
# create the single cell objects
        # the NearestCentroidClassifier works (not well).
        # We need the current markers in the foldTester, the log normalized
        # data to use with the classifiers, the rank data to create the next
data = Rocks(X, y.astype(int), verbose=1)
data.normalize(log=False, totalexpr=10000)
ft = FoldTester(data)

# currently assuming that we have some folds to load
ft.folds = folds

ft.makerocks(0)

# it's now safe to log-normalize the original rocks object since 
# all of the rocks objects in the ft made copies of the original data
# this will make classification work.
data.normalize(log=True, totalexpr=10000)

print("Data uploaded sucessfully!", flush=True)


########### Run some tests   ####################
# Data are saved to named files in the current working directory: check the script for more specifics
# parameters are now read in from files.  Uncomment these if you want to force things.

# small mesh
#lambMesh = list(np.linspace(0.01,0.1,num=10,endpoint=True))
#alphaMesh = (np.linspace(0.01, 0.1, num=10, endpoint=True))
# large mesh
#lambMesh = list(np.linspace(0.1,1.0,num=10,endpoint=True))
#alphaMesh = (np.linspace(0.05, 0.5, num=10, endpoint=True))

# store lambFac, alpha, #correct for each value of currLamb
ovrResults = np.zeros([len(lambs), len(alphaMesh)*len(lambMesh),3])
# store lists of number of support genes.  5 is the number of folds here.
numGenes = np.zeros([len(lambs), len(alphaMesh)*len(lambMesh), 5], dtype=int)

for eye in range(len(lambMesh)):
    
    lambFac = lambMesh[eye]
    
    for jay in range(len(alphaMesh)):
        
        alpha = alphaMesh[jay]
        
        # index in ovrResults
        index = len(alphaMesh)*eye + jay

        for kay in range(len(lambs)):

            currLamb = lambs[kay]
            print("Working on lamb = {}, alpha = {}, lambFac = {}".format(currLamb, alpha, lambFac))

            # This finds the data for our current parameters
            ft.selectmarkers(lambda sc: Rocks.findCSmarkers(sc, currLamb=currLamb, 
                alpha=alpha, lambFac=lambFac, writeOut = False), verbose = 1)

            # the actual classification
            ft.classify(RandomForest)
            currCorrect = y[ ft.yhat == y.flatten() ].shape[0]

            ovrResults[kay, index] = np.array([lambFac, alpha, currCorrect])
            numGenes[kay, index] = np.array([len(a) for a in ft.markers], dtype=int)
            np.savez("yhat-lamb{}-lfac{}-a{}.npz".format(currLamb, lambFac, alpha), yhat=ft.yhat)
    
# write the results to a file
fName = 'ovrResults-{}.dat'.format(arrID)
rFile = open(fName, 'wb')

fName = 'ovrNumGenes-{}.dat'.format(arrID)
nFile = open(fName, 'wb')


for kay in range(len(lambs)):
    comment = "Current lambda: " + str(lambs[kay])
    np.savetxt(rFile, ovrResults[kay], fmt=['%1.2f', '%1.2f', '%5i'], header=comment)
    np.savetxt(nFile, numGenes[kay], header=comment)

rFile.close()
nFile.close()
