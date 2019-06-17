04 June 2019

Files in this directory:

* `louvain.csv` contains a Louvain clustering of the 10x mouse data set that was generated as part of a scanpy tutorial.  See `https://github.com/theislab/scanpy_usage/tree/master/170522_visualizing_one_million_cells` commit `ba6eb85`.

* `10x-proc.ipynb` contains the initial processing of the 10x data.  We load the data in from the download (on the 10x website, see the paper for the link), add the clustering information, and delete all genes that are not expressed.  We also explore the data set and why it requires so much memory (~32GB, even in sparse format) to load.  We run the RankCorr method (known as 1bcs) here and save the vectors of dot products that we generate.  These files are called, for example, `1M-fold0-consts.npz`.  They are also in this directory.

The method that we used to parallelize the RankCorr method is very buggy, at least when used in a jupyter notebook.  Thus, we implement the scanpy methods in a python script and manually parallelize (by running each fold separately).  


* `10x-1bcs-markerProc.ipynb` takes the vectors of dot products that are created in `10x-proc.ipynb` and finishes the processing to get lists of markers. We save vectors of values of $s$ that would be needed to select each marker in `rc-svals-map.npz`.   (This is the final version that we use).


* `10x-markerProc.ipynb` processes all of the markers (from both scanpy methods and RankCorr) to create predictions for all the cells in various `yhat` files.  These files are not included here, as they are quite large.  Please contact the authors for this data if you would like them.

* `10x-graphs.ipynb` takes the yhat information and comapares to the original Louvain clustering to create the graphs that appear in the paper.  Without the yhat data, it will be difficult (or impossible) to run this.

* The directory `scanpy` contains the  file that was used to generate all of the data from the scanpy methods (Wilcoxon, t-test, and logisitc regression) as well as the data from these three methods in the fold directories.
