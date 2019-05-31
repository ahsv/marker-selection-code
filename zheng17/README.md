31 May 2019

* `zheng17-proc.ipynb` Contains the code used to process the Zheng data sets.  Starting with the data downloaded from the 10x website, it creates `annData` objects for both the ZhengFull and ZhengFilt data sets.  It also creats the folds needed for cross validation.  It also contains Louvain clustering and an exploration of sparsity.

* `zheng17-5folds.npz` contains the folds that were used for cross-validation on the ZhengFull and ZhengFilt data sets using the bulk labels.
