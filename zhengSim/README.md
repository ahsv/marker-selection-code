21 May 2019

The flow that was used to create the ZhengSim data.

1) `two-group-processing.ipynb` was run to create subsets of the b cell data set from Zheng 17 (the sample that was used to define the b Cell cluster in the bulk labels).  We choose to work with one sample since it matches quite well to the original 68k PBMC data set, while the second sample (of t cells) does not.  Thus, we do not combine the two samples and instead just use Splatter to simulate data based on the one sample.

2) The 20 data sets were genereated using Splatter in R 3.5 .  For a single data set, we run `splatSim.R` and then `procSplat.sh` to turn the results into data that can be easily loaded into python.  In particular, `splatSim.R` outputs `cellInfo.dat` which is turned into `y.dat` (the vector of cluster labels), `geneInfo.dat` which is turned into `deGenes.dat` (the genes that are differentially expressed), and `counts.dat` (that contains the counts).We include here the parameters that were learned (by splatter) for each data set that was generated; learning these parameters is the most computationally intensive portion of the simulation.

3) `splat-simulated-data.ipynb` is used to filter the genes from the `allGenes` simulation condition to create the "filter after simulation" data sets.  This information is saved in `filtered-splat.h5ad` and `filtered-splat-deInfo.npz` (the differentially expressed genes that pass through the filter).  We also sparsify the `allGenes` simulation conditions to make the files easier to work with (`sparseCounts.npz`).  Finally, the folds are different for each trial - thus, a set of folds is generated for each of the 20 data sets.

The simulated data is then processed in the following way:

1) In `1bcs.ipynb`, we apply the RankCorr algorithm to the simulated data sets.  We generate the vectors of dot products that we use for fast marker selection (files with `consts` in their name).  Then we also generate a cross-validated classification error (in an `errorRate` file) and non-cross validation precision, tpr, and fpr information (in a `statInfo` file). `1bcs.ipynb` is in the root of this repo.

2) In `scanpy-pvals.ipynb` (in the `scanpy` directory at the root of this repo), we generate cross-validated marker lists for each of the 30 data sets (stored in `markerList` files) and non-cross validated marker lists (for the generation of precision, tpr, and fpr information; stored in `valsNmarks` files).  Marker lists are generated for the Wilcoxon, t-test, and Logisitic Regression methods.

3) In `zhengsim-pval-classification-errors.ipynb` we process the wilcoxon, t-test, and logisitc regression markers to create corresponding cross validated classification error rate data (stored in `errorRate` files) out of the `markerList` files and `statInfo` information out of the `valsNmarks` files.  We also generate `ave-errorRate` files for each of the simulation conditions: these files contain the average classification error rates across all 10 trials for all four of the methods in consideration here.

4) In `zhengsim-graphs.ipynb`, we use the `ave-errorRate` files and the `statInfo` files to make the plot of classification error rate, precision, tpr, and fpr that appear in the paper.


We do not include all of the data here, but it is available on request.  
