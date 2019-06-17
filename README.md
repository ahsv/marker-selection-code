June 08 2019

This repo contains the code that is needed to reproduce the data presented in the manuscript "Comparision of marker selection methods for high thoughput scRNA-seq data."

The data sets themselves are not included here: see the manuscript for links to the data set.  You will need to change the paths (in the files) to load the correct data sets, where necessary.

* `loadData.py` contains methods for loading the data sets in multiple different formats.  Changing the paths in `loadData.py` should help a lot.

* The `classification` directory contains the notebooks that were used to compute the classification metrics and create the classification plots that are found in the manuscript.  

* The `clustering` directory contains the notebooks that were used to compute the clustering metrics and create the clustering plots that are found in the manuscript.  It also contains the summarized clustering data (that was generated from lists of markers).

* The directories named after data sets (`paul`, `zeisel`, `zheng17`, `10xmouse`, and `zhengSim`) contain the scripts that were used to initially process the data (if any) as well as any specific method implementations that were used to find markers on those data sets.  Those directories also (mostly) contain the markers that were selected by the different methods.

* The `scvi` and `enets` directories contain the implementations of and data collected using the scvi and elastic nets methods, respectively.  

* The `scanpy` directory contains the file that was used to find markers using the Wilcoxon, t-test, and logistic regression methods (these were implemented in the `scanpy` package).  It also contains the modified files that we used to fix a couple of bugs in scanpy.  See the paper for more details.

* The `picturedrocks` directory contains a modified old version of the picturedrocks class by Umang Varma.  These files also contain the implementation of Spa and RankCorr that were used for our paper. 
