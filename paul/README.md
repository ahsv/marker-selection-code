04 June 2019

### All processing of the Paul and Zeisel marker data is found in `PaulZeisel-classification-graphs.ipynb` (for the classification metrics) and `clustering-metrics.ipynb` (for the clustering metrics).

This directory mainly contains implementations of the methods as well as some raw marker data.

The files used to generate the classification pictures are `plotInfo` files.  They are generated using the marker data in the file `paulZeisel-classification-graphs.ipynb` which is found in the `classification` directory at the root of this repo.  The `plotInfo` files are found in the `plot` subdirectory.

* `d3e.py` contains the code needed to run d3e.  It should be put in the same directory as the d3e files (which might need to be converted to python 3, depending on the system that you are running on).  Also in the directory should be `paul-base.dat`: the paul data set, with each row containing the information from a gene, and one blank row at the top.  The data is read into the program through scanpy, so we need access to that.   Contact me if the supporting files are desired and see the paper for the computational resources required.

* for scVI, see the `scvi` directory located at the root of this repo.  The markers are found in that directory.
* for Elastic Nets, see the `enets` directory located at the root of this repo. The markers are found in that directory.
* For Wilcoxon, the t-test, and Logisitic Regression (all implemented in `scanpy`) see the `scanpy-methods` directory at the root of this repo.  The marker lists from these three methods are included here in the `markerList` files.
* For RankCorr, see the file `1bcs.ipynb` that is in the root of this repo.  It outputs lists of yhats (in `plotInfo`) files - we do not save any markers.
* Random markers are selected in the `Random_markers.ipynb` notebook in the root directory of this repository.


### The `edgeR` directory

Contains the code used to run edgeR on Paul.

* `edgeR.py` contains the edgeR script.  It is similar to all of the other edgeR scripts, specialized for the paul data set.
* `runTests.pbs` contains the file that was used to submit the job.  It has info about which packages were used as well as the resources that were requested.  

For an array job, you are expected to set up five directories (called `fold0`, ..., `fold4`).  In each directory, there should be a file called `foldNum.dat` that contains the fold number (alone, on one line).  The data for each fold will be collected in the specified directory if you are running on a torque pbs system.

The marker lists (p-value lists, processed in the clustering and classification metric files) are included here as well.

### The `mast` directory

Contains the code used to run  MAST on Paul.

* `mast.py` contains the MAST script.  It is similar to all of the other MAST scripts, specialized for the Paul data set.
* `runTests.pbs` contains the file that was used to submit the job.  It has info about which packages were used as well as the resources that were requested.  

For an array job, you are expected to set up five directories (called `fold0`, ..., `fold4`).  In each directory, there should be a file called `foldNum.dat` that contains the fold number (alone, on one line).  The data for each fold will be collected in the specified directory if you are running on a torque pbs system.

The marker lists (p-value lists, processed in the clustering and classification metric files) are included here as well.

### The `edgeRdet` and `mastDet` directories

Are similar to the above, for the methods that incorporate the detection rate.

### The `spa` directory 

Includes the files needed to run SPA.  For faster processing, we split the input mesh of paramters into 10 groups and ran 10 parallel jobs on those groups.  Each group was run in a separate directory and then the data was collected, processed, and compiled (using the `process_script.sh` file).  The scripts that we used are all in the `spa` directory, as well as an example input directory set up with a parameter input file.  We changed alpha from 0.01 to 0.1 in our different processing directories.

minMarks.dat contains the input parameters that resulted in the minimum number of markers selected for a given input value of s.

maxCorr.dat contains the input parameters that resulted in the maximum number of cells classified correctly according to the NCC for a given input value of s.

effMarks.dat contains the input parameters that resulted in the maximum ratio (number of cells correctly classified by the NCC)/(number of markers selected).  

maxCorr.dat tends to produce the best results, and thus those are the values that are reported.  The corresponding classification vectors (of cells) are included in `paul15-genzel-plotInfo.npz` and `paul15-genzel-plotInfo-RF.npz`; see the notebook `1bcs.ipynb` for more information about how these files were created.  The markers are selected for these parameters on-the-fly in the `clustering-metric.ipynb` processing notebook and the `plotInfo` files are in the `plot` directory..


