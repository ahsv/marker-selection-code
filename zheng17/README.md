31 May 2019

* `zheng17-proc.ipynb` Contains the code used to process the Zheng data sets.  Starting with the data downloaded from the 10x website, it creates `annData` objects for both the ZhengFull and ZhengFilt data sets.  It also creats the folds needed for cross validation.  It also contains Louvain clustering and an exploration of sparsity.

* `zheng17-5folds.npz` contains the folds that were used for cross-validation on the ZhengFull and ZhengFilt data sets using the bulk labels.

The `markerList` files contain the selected markers (in order) by the methods in the file names.  They were generated using the `scanpy-pvalues.ipynb` notebook: see the `scanpy` directory at the root of this repo.



### The `rankcorr` directory

* `1bcs-zheng.ipynb` contains the code needed to run RankCorr on the Zheng data sets.  Currently it pulls in the Rocks class but only uses it for the `sparse_dot_tau` function as far as I can tell.

The other files are created by `1bcs-zheng.ipynb` and include the lists of markers and lists of yhats that are used for the creation of the graphs in the paper.


### The `edgeR` directory

Contains the code used to run edgeR on ZhengFilt.

* `edgeR.py` contains the edgeR script.  It is similar to all of the other edgeR scripts, specialized for the zheng data set.
* `runTests.pbs` contains the file that was used to submit the job.  It has info about which packages were used as well as the resources that were requested.  

For an array job, you are expected to set up five directories (called `fold0`, ..., `fold4`).  In each directory, there should be a file called `foldNum.dat` that contains the fold number (alone, on one line).  The data for each fold will be collected in the specified directory if you are running on a torque pbs system.

The marker lists (p-value lists, processed in the clustering and classification metric files) are included here as well.

### The `mast` directory

Contains the code used to run  MAST on ZhengFilt.

* `mast.py` contains the MAST script.  It is similar to all of the other MAST scripts, specialized for the zheng data set.
* `runTests.pbs` contains the file that was used to submit the job.  It has info about which packages were used as well as the resources that were requested.  

For an array job, you are expected to set up five directories (called `fold0`, ..., `fold4`).  In each directory, there should be a file called `foldNum.dat` that contains the fold number (alone, on one line).  The data for each fold will be collected in the specified directory if you are running on a torque pbs system.

The marker lists (p-value lists, processed in the clustering and classification metric files) are included here as well.
